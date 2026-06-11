#!/bin/sh
# download bioreason-pro GO function test set (cafa5 temporal holdout) into data/bioreason_go
set -e
here="$(dirname "$0")"
d="$here/../data/bioreason_go"
hf="$here/../.venv/bin/hf"; [ -x "$hf" ] || hf=hf   # prefer repo venv, else PATH
"$hf" download wanglab/bioreason-pro-test-data --repo-type dataset --local-dir "$d/test"
