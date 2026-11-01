#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# != 2 ]]; then echo 'Usage: scripts/export_viewer.sh WORLD_ROOT OUTPUT_DIR' >&2; exit 2; fi
world=$(realpath "$1")
out=$(realpath -m "$2")
if [[ -e "$out" ]]; then echo 'Output exists; choose a new directory.' >&2; exit 2; fi
export CARGO_TARGET_DIR=${CARGO_TARGET_DIR:-"$root/target"}
export TMPDIR=${MAIFORGE_TMP_DIR:-"$CARGO_TARGET_DIR/tmp"}
export CARGO_BUILD_JOBS=${CARGO_BUILD_JOBS:-4}
mkdir -p "$TMPDIR" "$(dirname "$out")"
cargo run --locked --manifest-path "$root/worldfoundry/Cargo.toml" -p pf-lab --example maiforge_viewer_export -- "$world" "$out"
