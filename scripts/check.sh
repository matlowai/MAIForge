#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
export CARGO_TARGET_DIR=${CARGO_TARGET_DIR:-"$root/target"}
export TMPDIR=${MAIFORGE_TMP_DIR:-"$CARGO_TARGET_DIR/tmp"}
export CARGO_BUILD_JOBS=${CARGO_BUILD_JOBS:-4}
export RUST_MIN_STACK=${RUST_MIN_STACK:-16777216}
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$TMPDIR"
python3 -B -m unittest discover -s assets/foundry/botany/picea_abies_v0 -p 'test_*.py'
python3 -B -m unittest discover -s assets/foundry/botany/picea_growth_v1 -p 'test_*.py'
cargo test --locked --manifest-path worldfoundry/Cargo.toml -p pf-grid --lib
cargo test --locked --manifest-path worldfoundry/Cargo.toml -p pf-surface --test marine_fines
export PE_EVIDENCE_DIR=${PE_EVIDENCE_DIR:-"$CARGO_TARGET_DIR/producer-evidence-$(date +%Y%m%dT%H%M%S)"}
cargo test --locked --manifest-path worldfoundry/Cargo.toml -p pf-lab --test producer_evidence
cargo test --locked --manifest-path worldfoundry/crates/pf-rockspec/Cargo.toml --lib
cargo test --locked --manifest-path worldfoundry/crates/pf-hydro-catchment/Cargo.toml
