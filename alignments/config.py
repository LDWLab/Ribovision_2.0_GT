cif_path_share = ""
pdb_path_share = ""
chainid = ""
pdbid = ""

# Local KrakenD gateway that performs all external API calls. Override via the
# KRAKEND_BASE_URL environment variable in production if needed.
KRAKEND_BASE_URL = "http://127.0.0.1:8080"

# Filesystem cache for the external-API proxy (2-day TTL enforced in extapi.py).
EXTAPI_CACHE_PATH = "/tmp/extapi_cache"

# Persistent cache for deterministic local R2DT runs (no TTL by default).
R2DT_CACHE_PATH = "/tmp/r2dt_cache"