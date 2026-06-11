uv pip install --python .venv -r requirements.txt#!/bin/sh
# download proteingym dms substitutions into data/proteingym
set -e
d="$(dirname "$0")/../data/proteingym"
mkdir -p "$d"
curl -L -o "$d/DMS_substitutions.csv" \
  https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/main/reference_files/DMS_substitutions.csv
curl -L -o "$d/subs.zip" \
  https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.1/DMS_ProteinGym_substitutions.zip
cd "$d" && unzip -q -o subs.zip && rm subs.zip
