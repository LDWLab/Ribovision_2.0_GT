#!/usr/bin/env python3
"""Merge RiboVision's /rv/* endpoints into an existing global KrakenD config.

Idempotent and non-destructive: enables the root-level plugin loader and
replaces only endpoints whose path starts with ``/rv/`` while preserving every
other service's configuration. Writes the merged JSON to stdout.

Usage:
    merge_config.py <global_config.json> <ribovision-endpoints.json> <plugin_dir>
"""
import json
import sys


def merge(global_path, endpoints_path, plugin_dir):
    with open(global_path) as f:
        cfg = json.load(f)
    with open(endpoints_path) as f:
        rv_endpoints = json.load(f)

    folder = plugin_dir if plugin_dir.endswith("/") else plugin_dir + "/"
    cfg["plugin"] = {"pattern": ".so", "folder": folder}

    others = [
        e
        for e in cfg.get("endpoints", [])
        if not str(e.get("endpoint", "")).startswith("/rv/")
    ]
    cfg["endpoints"] = others + rv_endpoints
    return cfg


def main():
    if len(sys.argv) != 4:
        sys.stderr.write(__doc__)
        return 2
    cfg = merge(sys.argv[1], sys.argv[2], sys.argv[3])
    json.dump(cfg, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
