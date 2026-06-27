"""Persistent cache for deterministic local R2DT runs.

R2DT produces a 2D layout that is a pure function of the input sequence, and the
final response additionally depends on the structure inputs (entity id, chain id
and the CIF/PDB file or pdbid) used to compute the base pairs. We therefore keep
two cache levels:

  - layout cache : key = sha256(sequence)            -> raw ``r2dt.py draw`` JSON
  - result cache : key = sha256(seq|entity|chain|id)  -> merged final response

Entries never expire by default (the mapping is constant); set ``R2DT_CACHE_TTL``
(seconds) to enforce a maximum age.
"""

import os
import hashlib

import alignments.config as config
from alignments.disk_cache import DiskCache

R2DT_CACHE_PATH = os.environ.get(
    "R2DT_CACHE_PATH",
    getattr(config, "R2DT_CACHE_PATH", "/tmp/r2dt_cache"),
)
_ttl_env = os.environ.get("R2DT_CACHE_TTL")
R2DT_CACHE_TTL = int(_ttl_env) if _ttl_env else None  # None => never expire

_cache = DiskCache(R2DT_CACHE_PATH, ttl_seconds=R2DT_CACHE_TTL)


def _sha_text(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _layout_key(sequence):
    return "layout:" + _sha_text(sequence)


def _result_key(sequence, entity_id, chain_id, struct_identity):
    composite = "|".join(
        [sequence, str(entity_id), str(chain_id), str(struct_identity)]
    )
    return "result:" + _sha_text(composite)


def get_layout_bytes(sequence):
    hit = _cache.get(_layout_key(sequence))
    if hit is None:
        return None
    _meta, body = hit
    return body


def set_layout_bytes(sequence, body):
    _cache.set(_layout_key(sequence), body)


def get_result(sequence, entity_id, chain_id, struct_identity):
    return _cache.get_json(_result_key(sequence, entity_id, chain_id, struct_identity))


def set_result(sequence, entity_id, chain_id, struct_identity, obj):
    _cache.set_json(
        _result_key(sequence, entity_id, chain_id, struct_identity), obj
    )
