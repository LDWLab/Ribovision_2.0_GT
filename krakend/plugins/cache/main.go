// KrakenD http-client plugin providing a global, forced-TTL in-memory cache.
//
// Unlike KrakenD CE's built-in qos/http-cache (which only honours upstream
// Cache-Control headers), this plugin caches every successful GET response for
// a fixed TTL (default 2 days) regardless of what the upstream sends. The cache
// is process-wide, so a single KrakenD instance acts as a shared cache for all
// local services whose backends opt in via "plugin/http-client".
//
// Backend config:
//   "extra_config": {
//     "plugin/http-client": {
//       "name": "krakend-cache-client",
//       "ttl": 172800,            // optional, seconds (default 172800 = 2 days)
//       "cache_methods": ["GET"]  // optional, defaults to ["GET"]
//     }
//   }
package main

import (
	"bytes"
	"context"
	"io"
	"net/http"
	"strconv"
	"strings"
	"sync"
	"time"
)

func main() {}

// ClientRegisterer is the symbol KrakenD looks up to register http clients.
var ClientRegisterer = registerer("krakend-cache-client")

type registerer string

func (r registerer) RegisterClients(f func(
	name string,
	handler func(context.Context, map[string]interface{}) (http.Handler, error),
)) {
	f(string(r), r.registerClients)
}

const defaultTTL = 172800 * time.Second

type cacheEntry struct {
	status  int
	header  http.Header
	body    []byte
	expires time.Time
}

type cacheStore struct {
	mu    sync.RWMutex
	items map[string]cacheEntry
}

func (c *cacheStore) get(key string) (cacheEntry, bool) {
	c.mu.RLock()
	e, ok := c.items[key]
	c.mu.RUnlock()
	if !ok || time.Now().After(e.expires) {
		return cacheEntry{}, false
	}
	return e, true
}

func (c *cacheStore) set(key string, e cacheEntry) {
	c.mu.Lock()
	c.items[key] = e
	c.mu.Unlock()
}

func (c *cacheStore) cleanup() {
	now := time.Now()
	c.mu.Lock()
	for k, e := range c.items {
		if now.After(e.expires) {
			delete(c.items, k)
		}
	}
	c.mu.Unlock()
}

// A single process-wide store shared across every backend that uses the plugin.
var (
	store       = &cacheStore{items: make(map[string]cacheEntry)}
	httpClient  = &http.Client{Timeout: 60 * time.Second}
	janitorOnce sync.Once
)

func (r registerer) registerClients(_ context.Context, extra map[string]interface{}) (http.Handler, error) {
	ttl := defaultTTL
	cacheMethods := map[string]bool{http.MethodGet: true}

	cfg := pluginConfig(extra, string(r))
	if v, ok := cfg["ttl"]; ok {
		ttl = parseTTL(v, ttl)
	}
	if v, ok := cfg["cache_methods"]; ok {
		if methods, ok := v.([]interface{}); ok && len(methods) > 0 {
			cacheMethods = map[string]bool{}
			for _, m := range methods {
				if s, ok := m.(string); ok {
					cacheMethods[strings.ToUpper(s)] = true
				}
			}
		}
	}

	janitorOnce.Do(func() {
		go func() {
			ticker := time.NewTicker(10 * time.Minute)
			for range ticker.C {
				store.cleanup()
			}
		}()
	})

	return http.HandlerFunc(func(w http.ResponseWriter, req *http.Request) {
		cacheable := cacheMethods[req.Method]
		key := req.Method + " " + req.URL.String()

		if cacheable {
			if e, ok := store.get(key); ok {
				writeResponse(w, e.header, e.status, e.body, "HIT")
				return
			}
		}

		var bodyReader io.Reader
		if req.Body != nil {
			b, _ := io.ReadAll(req.Body)
			bodyReader = bytes.NewReader(b)
		}
		upReq, err := http.NewRequestWithContext(req.Context(), req.Method, req.URL.String(), bodyReader)
		if err != nil {
			http.Error(w, err.Error(), http.StatusBadGateway)
			return
		}
		copyHeader(upReq.Header, req.Header)
		// Let Go transparently handle compression so cached bodies are decoded.
		upReq.Header.Del("Accept-Encoding")

		resp, err := httpClient.Do(upReq)
		if err != nil {
			http.Error(w, err.Error(), http.StatusBadGateway)
			return
		}
		defer resp.Body.Close()
		respBody, _ := io.ReadAll(resp.Body)

		if cacheable && resp.StatusCode >= 200 && resp.StatusCode < 300 {
			store.set(key, cacheEntry{
				status:  resp.StatusCode,
				header:  resp.Header.Clone(),
				body:    respBody,
				expires: time.Now().Add(ttl),
			})
		}

		writeResponse(w, resp.Header, resp.StatusCode, respBody, "MISS")
	}), nil
}

func writeResponse(w http.ResponseWriter, header http.Header, status int, body []byte, cacheState string) {
	copyHeader(w.Header(), header)
	w.Header().Set("X-Krakend-Cache", cacheState)
	w.WriteHeader(status)
	w.Write(body)
}

func copyHeader(dst, src http.Header) {
	for k, vv := range src {
		for _, v := range vv {
			dst.Add(k, v)
		}
	}
}

// pluginConfig extracts this plugin's options from the backend extra_config,
// tolerating the various nesting shapes KrakenD may pass.
func pluginConfig(extra map[string]interface{}, name string) map[string]interface{} {
	if v, ok := extra["plugin/http-client"]; ok {
		if m, ok := v.(map[string]interface{}); ok {
			return m
		}
	}
	if v, ok := extra[name]; ok {
		if m, ok := v.(map[string]interface{}); ok {
			return m
		}
	}
	return extra
}

func parseTTL(v interface{}, def time.Duration) time.Duration {
	switch t := v.(type) {
	case float64:
		if t > 0 {
			return time.Duration(t) * time.Second
		}
	case int:
		if t > 0 {
			return time.Duration(t) * time.Second
		}
	case string:
		if n, err := strconv.Atoi(t); err == nil && n > 0 {
			return time.Duration(n) * time.Second
		}
	}
	return def
}
