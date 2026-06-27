"""Same-origin proxy for all external API calls.

The browser only ever talks to these ``/extapi/`` endpoints. Each view forwards
to the shared global KrakenD gateway under the ``/rv/`` namespace. That single
KrakenD instance is the global cache for every local service: it performs the
real external call and caches each successful GET for a forced 2-day TTL (via
the ``krakend-cache-client`` plugin), so this layer stays thin and just adds:
  - a lightweight per-IP rate limit,
  - graceful fallbacks that preserve the response shapes the frontend expects.

An optional local filesystem cache (off by default) can be enabled with
``EXTAPI_LOCAL_CACHE=1`` for extra resilience across gateway restarts.
"""

import os
import time
import logging
import threading
from functools import wraps
from urllib.parse import urlencode, quote

import requests

from django.http import HttpResponse, JsonResponse
from django.views.decorators.csrf import csrf_exempt

import alignments.config as config
from alignments.disk_cache import DiskCache

logger = logging.getLogger("ribovision3-logger")

KRAKEND_BASE_URL = os.environ.get(
    "KRAKEND_BASE_URL",
    getattr(config, "KRAKEND_BASE_URL", "http://127.0.0.1:8080"),
).rstrip("/")

# All RiboVision endpoints live under this namespace on the shared gateway.
RV_PREFIX = os.environ.get("KRAKEND_RV_PREFIX", "/rv").rstrip("/")

# The global KrakenD cache is authoritative. The local disk cache is opt-in.
LOCAL_CACHE_ENABLED = os.environ.get("EXTAPI_LOCAL_CACHE", "0") == "1"
EXTAPI_CACHE_TTL = int(os.environ.get("EXTAPI_CACHE_TTL", 172800))  # 2 days
EXTAPI_CACHE_PATH = os.environ.get(
    "EXTAPI_CACHE_PATH",
    getattr(config, "EXTAPI_CACHE_PATH", "/tmp/extapi_cache"),
)
DEFAULT_TIMEOUT = int(os.environ.get("EXTAPI_TIMEOUT", 60))

_cache = DiskCache(EXTAPI_CACHE_PATH, ttl_seconds=EXTAPI_CACHE_TTL)

# ---------------------------------------------------------------------------
# Per-IP rate limiting (lightweight in-memory sliding window).
# ---------------------------------------------------------------------------
RATE_LIMIT_MAX = int(os.environ.get("EXTAPI_RATE_LIMIT_MAX", 120))
RATE_LIMIT_WINDOW = int(os.environ.get("EXTAPI_RATE_LIMIT_WINDOW", 60))  # seconds
_rate_lock = threading.Lock()
_rate_state = {}


def _client_ip(request):
    forwarded = request.META.get("HTTP_X_FORWARDED_FOR")
    if forwarded:
        return forwarded.split(",")[0].strip()
    return request.META.get("REMOTE_ADDR", "unknown")


def rate_limited(view):
    @wraps(view)
    def wrapper(request, *args, **kwargs):
        ip = _client_ip(request)
        now = time.time()
        with _rate_lock:
            hits = [t for t in _rate_state.get(ip, []) if now - t < RATE_LIMIT_WINDOW]
            if len(hits) >= RATE_LIMIT_MAX:
                _rate_state[ip] = hits
                logger.warning("extapi rate limit exceeded for %s", ip)
                resp = JsonResponse(
                    {"error": "Too many requests"}, status=429
                )
                resp["Retry-After"] = str(RATE_LIMIT_WINDOW)
                return resp
            hits.append(now)
            _rate_state[ip] = hits
        return view(request, *args, **kwargs)

    return wrapper


# ---------------------------------------------------------------------------
# Core proxy helper.
# ---------------------------------------------------------------------------
def _build_url(krakend_path, query=None):
    url = KRAKEND_BASE_URL + RV_PREFIX + krakend_path
    if query:
        url += "?" + query
    return url


def _proxy(request, krakend_path, query=None, cacheable=True, timeout=DEFAULT_TIMEOUT):
    url = _build_url(krakend_path, query)
    cache_key = "%s %s" % (request.method, url)

    is_get = request.method == "GET"
    if cacheable and is_get and LOCAL_CACHE_ENABLED:
        hit = _cache.get(cache_key)
        if hit is not None:
            meta, body = hit
            resp = HttpResponse(
                body,
                status=meta.get("status", 200),
                content_type=meta.get("content_type", "application/octet-stream"),
            )
            resp["X-Extapi-Cache"] = "HIT"
            return resp

    try:
        if is_get:
            r = requests.get(url, timeout=timeout)
        else:
            content_type = request.META.get("CONTENT_TYPE", "")
            headers = {"Content-Type": content_type} if content_type else {}
            r = requests.request(
                request.method,
                url,
                data=request.body,
                headers=headers,
                timeout=timeout,
            )
    except requests.RequestException as exc:
        logger.error("extapi proxy failed for %s: %s", url, exc)
        return JsonResponse(
            {"error": "Upstream request failed", "detail": str(exc)},
            status=502,
        )

    content_type = r.headers.get("Content-Type", "application/octet-stream")
    body = r.content

    if cacheable and is_get and r.status_code == 200 and LOCAL_CACHE_ENABLED:
        try:
            _cache.set(
                cache_key,
                body,
                {"status": r.status_code, "content_type": content_type},
            )
        except OSError as exc:
            logger.warning("extapi cache write failed for %s: %s", url, exc)

    resp = HttpResponse(body, status=r.status_code, content_type=content_type)
    resp["X-Extapi-Cache"] = "MISS"
    return resp


# ---------------------------------------------------------------------------
# Endpoint views.
# ---------------------------------------------------------------------------
@rate_limited
def pdbe_molecules(request, pdb):
    return _proxy(request, "/pdbe/molecules/%s" % quote(pdb.lower()))


@rate_limited
def pdbe_polymer_coverage(request, pdb, chain):
    return _proxy(
        request,
        "/pdbe/polymer-coverage/%s/%s" % (quote(pdb), quote(chain)),
    )


@rate_limited
def pdbe_static_entry(request, pdb, entity, chain):
    return _proxy(
        request,
        "/pdbe/static-entry/%s/%s/%s"
        % (quote(pdb.lower()), quote(entity), quote(chain)),
    )


@rate_limited
def rcsb_model(request, pdb):
    query = request.META.get("QUERY_STRING", "")
    return _proxy(request, "/rcsb/model/%s" % quote(pdb.lower()), query=query)


@rate_limited
def pdbe_model_server(request, pdb):
    query = request.META.get("QUERY_STRING", "")
    return _proxy(request, "/pdbe/model-server/%s" % quote(pdb.lower()), query=query)


@rate_limited
def ribosome_banclass(request, pdb, chain):
    query = urlencode({"pdbid": pdb, "auth_asym_id": chain, "format": "json"})
    return _proxy(request, "/ribosome/banclass", query=query)


@rate_limited
def rcsb_graphql(request):
    query = request.META.get("QUERY_STRING", "")
    return _proxy(request, "/rcsb/graphql", query=query)


@rate_limited
def bgsu_basepairs(request, pdb, chain):
    return _proxy(
        request,
        "/bgsu/basepairs/%s/%s" % (quote(pdb.lower()), quote(chain)),
    )


@rate_limited
def bgsu_basepairs_nested(request, pdb, chain):
    return _proxy(
        request,
        "/bgsu/basepairs-nested/%s/%s" % (quote(pdb.lower()), quote(chain)),
    )


# BLAST is job-specific and intentionally not cached.
@csrf_exempt
@rate_limited
def blast_run(request):
    return _proxy(request, "/blast/run", cacheable=False)


@rate_limited
def blast_status(request, jobid):
    return _proxy(request, "/blast/status/%s" % quote(jobid), cacheable=False)


@rate_limited
def blast_result(request, jobid):
    return _proxy(request, "/blast/result/%s" % quote(jobid), cacheable=False)
