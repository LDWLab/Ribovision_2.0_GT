#!/bin/bash
# Builds the KrakenD forced-TTL cache plugin.
# Must be built with the SAME Go version KrakenD was built with (go1.20.14 for
# KrakenD 2.6.0) or KrakenD will refuse to load the .so.
set -e

GO="${GO:-/usr/local/go/bin/go}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${1:-$SCRIPT_DIR/cache.so}"

echo "Using $($GO version)"
cd "$SCRIPT_DIR"
CGO_ENABLED=1 "$GO" build -buildmode=plugin -o "$OUT" .
echo "Built plugin: $OUT"
