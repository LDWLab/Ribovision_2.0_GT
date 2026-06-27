"""Per-environment filesystem paths (FR3D / R2DT / logs).

These are environment-specific and live OUTSIDE version control. They are read
from the same sources as the Django settings, where later sources override
earlier ones:
  1. /etc/ribovision_config.json    (system-wide, e.g. the production server)
  2. DESIRE/local_config.json       (repo-local, gitignored)
  3. process environment variables
Copy DESIRE/local_config.json.example -> DESIRE/local_config.json and fill in
the FR3D_PATH / R2DT_PATH / LOGS_PATH values. See DEVELOPMENT.md.
"""
import os
import json

_BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load_env_config():
    merged = {}
    for path in ('/etc/ribovision_config.json',
                 os.path.join(_BASE_DIR, 'DESIRE', 'local_config.json')):
        try:
            with open(path) as fh:
                merged.update(json.load(fh))
        except (IOError, OSError, ValueError):
            pass
    return merged


_CONFIG = _load_env_config()


def _conf(key, default=''):
    value = _CONFIG.get(key)
    if value is None:
        value = os.environ.get(key, default)
    return value


FR3D_PATH = _conf('FR3D_PATH')
R2DT_PATH = _conf('R2DT_PATH')
LOGS_PATH = _conf('LOGS_PATH')