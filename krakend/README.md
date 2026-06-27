# Shared global KrakenD gateway

Every third-party API call (EBI PDBe, RCSB, ribosome.xyz, BGSU FR3D, EBI NCBI
BLAST) goes through **one shared, global KrakenD instance** on this host. The
browser never calls external services directly: it calls same-origin Django
`/extapi/...` endpoints (see `alignments/extapi.py`), Django forwards to the
gateway on `127.0.0.1:8080` under the `/rv/` namespace, and KrakenD performs the
real external call and caches it.

This gateway is **shared infrastructure** used by multiple local services
(e.g. the "Exornata API Gateway"). It is NOT owned by the RiboVision app deploy:

- Live config: `/etc/krakend/krakend.json`
- Systemd unit: `krakend` (`/usr/lib/systemd/system/krakend.service`)
- Plugins: `/etc/krakend/plugins/`

## Global forced-TTL cache (Go plugin)

`plugins/cache/` is a KrakenD `plugin/http-client` (`krakend-cache-client`) that
provides a **process-wide, forced 2-day in-memory cache**. Unlike KrakenD CE's
built-in `qos/http-cache` (which only honours upstream `Cache-Control`), this
plugin caches every successful GET for a fixed TTL regardless of upstream
headers. Because it is one process, the cache is shared by every service whose
backends opt in via `plugin/http-client`.

The plugin must be built with the **same Go version as KrakenD** (go1.20.14 for
KrakenD 2.6.0) or it will not load:

```bash
bash plugins/cache/build.sh        # produces plugins/cache/cache.so
krakend test-plugin -c plugins/cache/cache.so
```

Per-backend usage (already set on every `/rv/*` backend):

```json
"extra_config": {
  "plugin/http-client": { "name": "krakend-cache-client", "ttl": 172800 }
}
```

The cache is in-memory and is cleared on gateway restart. For an extra
persistent layer, Django can also cache locally by setting `EXTAPI_LOCAL_CACHE=1`
(off by default; the global gateway cache is authoritative).

## Install / update RiboVision on the global gateway

`install-ribovision-gateway.sh` is idempotent and non-destructive. It builds and
installs the plugin, then **merges** RiboVision's `/rv/*` endpoints into the
existing global config (via `merge_config.py`) — preserving every other
service's endpoints — backs up the live config, validates, and restarts:

```bash
bash krakend/install-ribovision-gateway.sh   # requires sudo
```

`ribovision-endpoints.json` is the source of RiboVision's endpoint definitions.
To add another service to the global cache, add its endpoints to
`/etc/krakend/krakend.json` (using the same `plugin/http-client` block) and
restart `krakend`.

## Endpoints (namespaced `/rv/`)
`/rv/pdbe/molecules/{pdb}`, `/rv/pdbe/polymer-coverage/{pdb}/{chain}`,
`/rv/pdbe/static-entry/{pdb}/{entity}/{chain}`, `/rv/pdbe/model-server/{pdb}`,
`/rv/rcsb/model/{pdb}`, `/rv/rcsb/graphql`, `/rv/ribosome/banclass`,
`/rv/bgsu/basepairs/{pdb}/{chain}`, `/rv/bgsu/basepairs-nested/{pdb}/{chain}`,
`/rv/blast/run`, `/rv/blast/status/{jobid}`, `/rv/blast/result/{jobid}`.

BLAST endpoints are intentionally uncached (job-specific).
