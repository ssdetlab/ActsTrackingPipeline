#!/usr/bin/env bash

# Absolute path to Run/ directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Build dir: three levels up from Run/, then build/ActsTrackingPipeline
# Run -> ActsTrackingPipeline -> src -> acts-work
BUILD_DIR="${SCRIPT_DIR}/../../../build/ActsTrackingPipeline"
BIN="${BUILD_DIR}/bin/RunPipeline"

# Config argument is relative to SCRIPT_DIR (Run/)
CONF_ARG="${1:-Confs/Pipeline.conf}"
CONF="${SCRIPT_DIR}/${CONF_ARG}"

if [ ! -x "$BIN" ]; then
  echo "Error: binary '$BIN' not found. Build with: cmake --build \"$BUILD_DIR\"" >&2
  exit 1
fi

if [ ! -f "$CONF" ]; then
  echo "Error: config '$CONF' not found" >&2
  exit 1
fi

# Run with cwd = Run/
cd "$SCRIPT_DIR" || exit 1
exec "$BIN" "$CONF"
