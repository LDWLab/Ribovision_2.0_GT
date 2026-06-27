#!/bin/bash

set -e

# Must match the directory Apache/mod_wsgi actually serves. The enabled vhost
# /etc/httpd/sites-available/ribovision3-ssl.conf points its python-path and
# WSGIScriptAlias at /home/github_repos/Ribovision_2.0_GT, so deploy here.
PROJECT_DIR=$PWD
OWNER="apache:apache"

cd "$PROJECT_DIR"

echo "=== Starting Deployment ==="

# 0. Cache sudo credentials upfront so later steps don't block
echo "[0/7] Requesting sudo credentials..."
sudo -v

# 1. Activate virtual environment
echo "[1/7] Activating Python environment..."
source env/bin/activate

# 2. Build react-msa-viewer
echo "[2/7] Building react-msa-viewer..."
cd react-msa-viewer/
NODE_ENV=production npx rollup -c
cd ..

# 3. Build pdbe-rna-viewer
echo "[3/7] Building pdbe-rna-viewer..."
cd pdbe-rna-viewer/
npm run build
cd ..

# 4. Build main project (webpack)
echo "[4/7] Building main project assets..."
npm run build

# 5. Collect static files
echo "[5/7] Collecting static files..."
python manage.py collectstatic --noinput

# 6. Fix ownership and permissions (BEFORE reloading the server)
# Apache runs as `apache`; every served file (incl. newly added Python modules)
# must be owned/readable by it, or imports fail silently and routes 404.
echo "[6/7] Fixing permissions..."
sudo chown -R "$OWNER" "$PROJECT_DIR"
sudo chmod -R o+r "$PROJECT_DIR"/static/
# Ensure all Python sources are group/owner readable for the apache worker.
sudo find "$PROJECT_DIR" -name '*.py' -exec chmod u+r,g+r {} +

# 7. Reload the server
echo "[7/7] Reloading WSGI..."
sudo touch DESIRE/wsgi.py

# NOTE: The external-API gateway is a SHARED, GLOBAL KrakenD instance (systemd
# unit `krakend`, config /etc/krakend/krakend.json) used by multiple services.
# It is intentionally NOT managed by this app deploy. To (re)install or update
# RiboVision's /rv/* endpoints + cache plugin on the global gateway, run:
#   bash krakend/install-ribovision-gateway.sh
echo "Note: global KrakenD gateway is managed separately (krakend/install-ribovision-gateway.sh)."

echo "=== Deployment Complete ==="