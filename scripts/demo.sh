#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
out=${1:-"$root/generated/demo-$(date +%Y%m%dT%H%M%S)"}
if [[ -e "$out" ]]; then echo 'Output exists; choose a new directory.' >&2; exit 2; fi
export CARGO_TARGET_DIR=${CARGO_TARGET_DIR:-"$root/target"}
export TMPDIR=${MAIFORGE_TMP_DIR:-"$CARGO_TARGET_DIR/tmp"}
export CARGO_BUILD_JOBS=${CARGO_BUILD_JOBS:-4}
mkdir -p "$TMPDIR"
cargo build --locked --manifest-path worldfoundry/Cargo.toml -p pf-lab -p pf-region -p pf-detail --bins
cargo build --locked --manifest-path worldfoundry/Cargo.toml -p pf-lab --example maiforge_inputs
"$CARGO_TARGET_DIR/debug/examples/maiforge_inputs" "$out"
mkdir -p "$out/logs"
"$CARGO_TARGET_DIR/debug/pf-region-skeleton" "$out/run" "$out/subgrid_relief_m.f32" "$out/site.json" "$out/t1a" > "$out/logs/skeleton.log" 2>&1
# Explicit bounded example: one 20 kyr local relaxation step, NOT equilibrium.
"$CARGO_TARGET_DIR/debug/pf-region-relief" run "$out/run" "$out/subgrid_relief_m.f32" "$out/site.json" "$out/t1a" "$out/t1b" simtime_kyr=20 dt_kyr=20 0,0 > "$out/logs/relief.log" 2>&1
python3 - "$out" <<'PY'
import json,sys
from pathlib import Path
r=Path(sys.argv[1]).resolve()
g=json.loads((r/'generation.json').read_text())
c={'t1b_window_dir':str(r/'t1b/w_0_0'),'t1a_dir':str(r/'t1a'),'t0_run_root':str(r),'site_record':str(r/'site.json'),'world_seed':g['seed'],'chart':0}
(r/'detail-config.json').write_text(json.dumps(c,indent=2)+'\n')
PY
"$CARGO_TARGET_DIR/debug/pf-detail" block "$out/detail-config.json" "$out/t2" 0 0 1 1 steps=8 > "$out/logs/detail.log" 2>&1
python3 scripts/export_terrain.py "$out"
echo "Generated example at $out. Read demo.json: file generation is not terrain-quality acceptance."
