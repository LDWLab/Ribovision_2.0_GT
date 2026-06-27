#!/bin/bash

set -e

PROJECT_DIR="/home/RiboVision3"
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
echo "[6/7] Fixing permissions..."
sudo chown -R "$OWNER" "$PROJECT_DIR"
sudo chmod -R o+r "$PROJECT_DIR"/static/

# 7. Reload the server
echo "[7/7] Reloading WSGI..."
sudo touch DESIRE/wsgi.py

echo "=== Deployment Complete ==="