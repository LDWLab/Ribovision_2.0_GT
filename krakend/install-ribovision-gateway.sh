#!/bin/bash
# Installs RiboVision into the EXISTING global KrakenD gateway.
#
# This is shared infrastructure: a single KrakenD instance (the "Exornata API
# Gateway" at /etc/krakend/krakend.json, systemd unit `krakend`) acts as the
# global cache for all local services. This script:
#   1. builds + installs the forced-TTL cache plugin into /etc/krakend/plugins/
#   2. enables the root-level "plugin" loader block (idempotent)
#   3. merges RiboVision's /rv/* endpoints into the global config (idempotent:
#      any existing /rv/* endpoints are replaced, every other service is kept)
#   4. validates the merged config and restarts the gateway
#
# It NEVER touches other services' endpoints. A timestamped backup of the global
# config is written before any change. Requires sudo.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GLOBAL_CONFIG="/etc/krakend/krakend.json"
PLUGIN_DIR="/etc/krakend/plugins"
PLUGIN_SO="$SCRIPT_DIR/plugins/cache/cache.so"
ENDPOINTS="$SCRIPT_DIR/ribovision-endpoints.json"
KRAKEND_USER="krakend"

echo "=== Installing RiboVision into the global KrakenD gateway ==="

# 1. Build the cache plugin if it is missing.
if [[ ! -f "$PLUGIN_SO" ]]; then
    echo "[1/6] Building cache plugin..."
    bash "$SCRIPT_DIR/plugins/cache/build.sh"
else
    echo "[1/6] Cache plugin already built: $PLUGIN_SO"
fi

# 2. Verify the plugin is loadable by this KrakenD build.
echo "[2/6] Verifying plugin loadability..."
krakend test-plugin -c "$PLUGIN_SO"

# 3. Install the plugin into the global plugin folder.
echo "[3/6] Installing plugin into $PLUGIN_DIR..."
sudo mkdir -p "$PLUGIN_DIR"
sudo cp "$PLUGIN_SO" "$PLUGIN_DIR/cache.so"
sudo chown -R "$KRAKEND_USER":"$KRAKEND_USER" "$PLUGIN_DIR" 2>/dev/null || true
sudo chmod 755 "$PLUGIN_DIR"
sudo chmod 644 "$PLUGIN_DIR/cache.so"

# 4. Back up + merge config (Python does the JSON surgery, never clobbers others).
echo "[4/6] Merging /rv/* endpoints into $GLOBAL_CONFIG..."
TS="$(date +%Y%m%d_%H%M%S)"
sudo cp "$GLOBAL_CONFIG" "${GLOBAL_CONFIG}.bak.${TS}"
echo "      backup: ${GLOBAL_CONFIG}.bak.${TS}"

MERGED="$(mktemp --suffix=.json)"
python3 "$SCRIPT_DIR/merge_config.py" "$GLOBAL_CONFIG" "$ENDPOINTS" "$PLUGIN_DIR" > "$MERGED"

# 5. Validate before going live.
echo "[5/6] Validating merged config..."
CHECK_TMP="$(mktemp -d)/krakend.json"
cp "$MERGED" "$CHECK_TMP"
# Plugin folder must be reachable during check; point it at the temp dir's sibling
# is unnecessary because the merged config references the real /etc plugin dir.
krakend check -c "$CHECK_TMP"

# 6. Install + restart.
echo "[6/6] Installing merged config and restarting gateway..."
sudo cp "$MERGED" "$GLOBAL_CONFIG"
sudo chown "$KRAKEND_USER":"$KRAKEND_USER" "$GLOBAL_CONFIG" 2>/dev/null || true
rm -f "$MERGED"
sudo systemctl restart krakend
sleep 2
sudo systemctl --no-pager --lines=5 status krakend || true

echo "=== Done. RiboVision endpoints are live under /rv/* on the global gateway. ==="
