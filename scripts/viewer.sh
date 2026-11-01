#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ ! -f "$root/viewer/package.json" ]]; then
  echo 'Initialize the adapter: git submodule update --init --recursive' >&2
  exit 2
fi
cd "$root/viewer"
npm ci
exec npm run dev -- "$@"
