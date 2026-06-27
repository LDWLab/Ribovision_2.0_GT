# Development & Deployment Guide

This repo uses **two long-lived branches**:

| Branch   | Purpose | External API calls | KrakenD / caching | Deploy configs |
|----------|---------|--------------------|-------------------|----------------|
| `master` | Local development. Runs as-is. | Direct (no gateway) | No | No |
| `proxy`  | Production server. | Routed through `/extapi` -> KrakenD gateway | Yes (gateway cache + backend caches) | Yes (`krakend/`, deploy perms) |

`proxy` is conceptually **`master` + a deployment layer** (KrakenD gateway, `/extapi`
routing, on-disk caches). All general feature work and bug fixes are made on
`master` and flow **`master -> proxy`**. KrakenD/caching code stays only on `proxy`.

---

## 1. Environment-specific config & secrets

**Nothing environment-specific is committed.** Secrets (Django `SECRET_KEY`, DB
credentials) and per-machine paths (`FR3D_PATH`, `R2DT_PATH`, `LOGS_PATH`) are read
at startup from these sources, where **later sources override earlier ones**:

1. `/etc/ribovision_config.json` — system-wide (used on the production server).
2. `DESIRE/local_config.json` — repo-local, **gitignored**.
3. Process environment variables.

The loader lives in `DESIRE/settings.py` (`_load_env_config`/`_conf`) and
`alignments/paths.py`. Recognized keys are documented in
`DESIRE/local_config.json.example`:

```
SECRET_KEY, DB_NAME, DB_ENGINE, DB_USER, DB_PASSWORD, DB_HOST, DB_PORT,
FR3D_PATH, R2DT_PATH, LOGS_PATH
```

Because `DESIRE/local_config.json` is gitignored, **each environment keeps its own
copy** and it is never overwritten by `git pull`/`merge`/`checkout`. That is how
dev and prod secrets/paths stay separate and preserved.

### First-time setup of a checkout
```bash
cp DESIRE/local_config.json.example DESIRE/local_config.json
# edit DESIRE/local_config.json with this environment's values
```
On the production server, `SECRET_KEY` already comes from
`/etc/ribovision_config.json`, so `DESIRE/local_config.json` there only needs the
DB credentials and paths.

---

## 2. Local development (branch `master`)

```bash
git clone <repo> && cd Ribovision_2.0_GT
git checkout master

# Python env
python3 -m venv env
source env/bin/activate
pip install -r requirements.txt          # includes mmcif-pdbx (needed by FR3D for CIF)

# Config
cp DESIRE/local_config.json.example DESIRE/local_config.json
# fill in SECRET_KEY, DB_*, FR3D_PATH, R2DT_PATH, LOGS_PATH

# Frontend bundles (needed for the viewer fixes that live in compiled JS)
npm install
npm run build                            # main project assets
( cd pdbe-rna-viewer && npm install && npm run build )
( cd react-msa-viewer && npm install && npm run build )

# Run
python manage.py collectstatic --noinput
python manage.py runserver
```

Quick sanity check that config loads and the DB is reachable:
```bash
python -c "import os,django;os.environ.setdefault('DJANGO_SETTINGS_MODULE','DESIRE.settings');django.setup();from django.db import connection;connection.ensure_connection();print('DB OK')"
```

### Notes
- `master` calls external APIs **directly**; no KrakenD gateway is required locally.
- Source fixes in Python and in directly-served JS (`static/alignments/RV3_helpers.js`,
  `static/js/components/*.js`) take effect immediately. Fixes in `*.vue` and
  `pdbe-rna-viewer/src/*.ts` require a build (`npm run build`) to appear in the
  served bundles.

---

## 3. Promoting changes to production (`master -> proxy`)

Do this on the production server (or any checkout), and never run it against the
live working tree mid-request without expecting a brief reload.

```bash
git checkout master && git pull
git checkout proxy  && git pull
git merge master                 # bring dev changes into production
```

**Expected merge conflict:** `static/js/components/getStructMappingAndTWC.js`
— `master` uses `ajax(...)`, `proxy` uses the cached `ajaxMappingCached(...)`.
Resolve by **keeping the `proxy` (cached) version** and applying any logic change
from `master` around it.

Then deploy:
```bash
bash deploy.sh                   # builds assets, collects static, reloads WSGI
```

KrakenD gateway changes are deployed separately and only when needed:
```bash
bash krakend/install-ribovision-gateway.sh
```

### Hotfix made directly on prod?
If you ever fix something directly on `proxy`, port it back to `master` with
`git cherry-pick <sha>` (from a `master` checkout/worktree) so the two branches
do not drift. Use a `git worktree` to avoid disturbing the live `proxy` checkout:
```bash
git worktree add ../rv-master-wt master
cd ../rv-master-wt && git cherry-pick <sha> && git push origin master
cd - && git worktree remove ../rv-master-wt
```

---

## 4. What lives where (quick reference)

- **Only on `proxy`:** `alignments/extapi.py`, `alignments/r2dt_cache.py`,
  `krakend/`, the `/extapi/*` routes in `alignments/urls.py`, the cached
  `ajax*Cached` helpers, deploy permission steps in `deploy.sh`.
- **On both branches:** all application logic, models, viewers, bug fixes,
  `DESIRE/settings.py` + `alignments/paths.py` (secret-free loaders),
  `DESIRE/local_config.json.example`.
- **Never committed (per environment):** `DESIRE/local_config.json`,
  `alignments/static/...` build caches, `.cache/`, `env/`.
