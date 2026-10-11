#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 scripts/check_torus_raw.py
python3 scripts/check_torus.py
python3 scripts/check_constants.py
python3 scripts/check_k3.py
python3 scripts/check_lean.py
