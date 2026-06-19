#!/bin/sh
# download BindingDB (curated protein-ligand binding affinities, articles subset) into data/bindingdb.
# the tsv's target sequence column holds each protein; affinity columns hold Ki/Kd/IC50/EC50.
set -e
d="$(dirname "$0")/bindingdb"
mkdir -p "$d"
curl -fsSL https://www.bindingdb.org/rwd/bind/downloads/BindingDB_BindingDB_Articles_202606_tsv.zip -o "$d/bindingdb.zip"
cd "$d" && unzip -q -o bindingdb.zip && rm bindingdb.zip
