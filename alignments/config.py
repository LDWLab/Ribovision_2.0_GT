import os as _os

cif_path_share = ""
pdb_path_share = ""
chainid = ""
pdbid = ""

# Local KrakenD gateway that performs all external API calls. Override via the
# KRAKEND_BASE_URL environment variable in production if needed.
KRAKEND_BASE_URL = "http://127.0.0.1:8080"

# Persistent cache root. NOTE: do NOT use /tmp here -- the Apache/httpd service
# runs with systemd PrivateTmp=true, so /tmp is a per-boot private mount that is
# WIPED on every httpd restart. We instead keep caches inside the project (one
# cache per environment: staging vs production never share) so deterministic
# results (R2DT layouts, MAFFT mappings) survive restarts. Override the root
# with the RV_CACHE_ROOT environment variable if a different location is wanted.
_PROJECT_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_CACHE_ROOT = _os.environ.get("RV_CACHE_ROOT", _os.path.join(_PROJECT_ROOT, ".cache"))

# Filesystem cache for the external-API proxy (2-day TTL enforced in extapi.py).
EXTAPI_CACHE_PATH = _os.path.join(_CACHE_ROOT, "extapi_cache")

# Persistent cache for deterministic local R2DT runs.
R2DT_CACHE_PATH = _os.path.join(_CACHE_ROOT, "r2dt_cache")
R2DT_CACHE_TTL = 604800  # 7 days

# Replace R2DT's canonical base pairs with FR3D-derived ones in the 2D viewer.
#
# FR3D indexes its base pairs by CIF label_seq_id (verified: 100% of endpoint
# nucleotide letters match _pdbx_poly_seq_scheme.mon_id at the same seq_id).
# This requires parse_cif4 to emit RNA_2D_json's `label_seq_ids` aligned to the
# R2DT drawing order, since the viewer keys nucleotide positions off
# label_seq_ids[i - 1] for svg_paths[i]. parse_cif4 does that alignment now; it
# previously zipped the drawing against _pdbx_poly_seq_scheme by ordinal
# position, which drifted whenever R2DT drew only part of the entity and put
# every base pair on the wrong nucleotide.
#
# Set FR3D_REPLACE_BASE_PAIRS=0 to fall back to R2DT's own canonical pairs.
FR3D_REPLACE_BASE_PAIRS = _os.environ.get("FR3D_REPLACE_BASE_PAIRS", "1") == "1"

# Persistent cache for deterministic MAFFT structure<->alignment mappings
# (/mapSeqAln/ and /mapSeqAlnOrig/). These outputs depend only on their inputs,
# so they never expire (TTL = None).
MAPPING_CACHE_PATH = _os.path.join(_CACHE_ROOT, "mapping_cache")

# Cache for DB-derived alignment builds (/ortholog-aln-api, /paralog-aln-api).
# These depend on the DESIRE database, so a TTL is used and the cache should be
# flushed when the database is repopulated (delete the directory).
ALN_CACHE_PATH = _os.path.join(_CACHE_ROOT, "aln_cache")
ALN_CACHE_TTL = 172800  # 2 days

# Persistent cache for raw structure coordinates fetched from external model
# servers (models.rcsb.org / litemol) by /custom-struc-data/. These coordinates
# are immutable for a given (pdb, entity) URL, so they never expire (TTL = None).
# This removes the "Loading Structure Data" network wait for repeat visitors and
# survives both Django session expiry and httpd restarts.
STRUCT_CACHE_PATH = _os.path.join(_CACHE_ROOT, "struct_cache")