"""Small filesystem cache used by the external-API proxy and the R2DT cache.

Stores each entry as two sharded files under ``base_dir``:
  - ``<hash>.meta`` : small JSON blob with at least a ``ts`` timestamp
  - ``<hash>.body`` : the raw response bytes

Writes are atomic (temp file + ``os.replace``) so concurrent readers never see a
partially written entry. A ``ttl_seconds`` of ``None`` means entries never expire
(used for the deterministic R2DT cache); a numeric value enforces a max age.
"""

import os
import json
import time
import hashlib
import tempfile


def _hash_key(key):
    return hashlib.sha256(key.encode("utf-8")).hexdigest()


def _atomic_write(path, data):
    directory = os.path.dirname(path)
    os.makedirs(directory, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(dir=directory)
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(data)
        os.replace(tmp_path, path)
    except Exception:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        raise


class DiskCache(object):
    def __init__(self, base_dir, ttl_seconds=None):
        self.base_dir = base_dir
        self.ttl_seconds = ttl_seconds
        try:
            os.makedirs(self.base_dir, exist_ok=True)
        except OSError:
            pass

    def _paths(self, key):
        h = _hash_key(key)
        shard = os.path.join(self.base_dir, h[:2])
        return (
            os.path.join(shard, h + ".meta"),
            os.path.join(shard, h + ".body"),
        )

    def get(self, key):
        """Return ``(meta_dict, body_bytes)`` on a fresh hit, else ``None``."""
        meta_path, body_path = self._paths(key)
        try:
            with open(meta_path, "r") as f:
                meta = json.load(f)
        except (OSError, ValueError):
            return None

        if self.ttl_seconds is not None:
            if (time.time() - float(meta.get("ts", 0))) > self.ttl_seconds:
                return None

        try:
            with open(body_path, "rb") as f:
                body = f.read()
        except OSError:
            return None

        return meta, body

    def set(self, key, body, meta=None):
        meta = dict(meta or {})
        meta["ts"] = time.time()
        meta_path, body_path = self._paths(key)
        # Write body first so a present meta always has a complete body.
        _atomic_write(body_path, body)
        _atomic_write(meta_path, json.dumps(meta).encode("utf-8"))

    def get_json(self, key):
        hit = self.get(key)
        if hit is None:
            return None
        _meta, body = hit
        try:
            return json.loads(body.decode("utf-8"))
        except (ValueError, UnicodeDecodeError):
            return None

    def set_json(self, key, obj, meta=None):
        self.set(key, json.dumps(obj).encode("utf-8"), meta)
