# protein-alignment-benchmarks

leakage-free qfit sequences to align a sequence model on, plus three benchmarks
to test if alignment improves prediction: thermompnn ddG, proteingym dms fitness,
and bioreason-pro GO function (cafa5 setup).

leakage guarantee: no sequence in an `align*.csv` is >20% identical to any
sequence in the **test set** of the benchmark it's "safe" from, so aligning on
it can't leak that benchmark's eval. we only filter against test seqs (eval is
test-only); leakage into a benchmark's train split is irrelevant and filtering
it just over-prunes qfit. one safe set per benchmark, plus one safe from all
three (kept / 37,038 qfit seqs):

- `align_ddg_safe.csv` — safe from the ddG test set (35,466)
- `align_dms_safe.csv` — safe from the dms set (24,490; proteingym is eval-only)
- `align_go_safe.csv` — safe from the GO test set (18,010)
- `align.csv` — safe from all three test sets (13,667)

## setup

```bash
uv venv --python 3.11 .venv
uv pip install --python .venv -r requirements.txt
sh src/download_ddg.sh                                               # ddG data
sh src/download_proteingym.sh                                        # dms data (~1GB)
sh src/download_bioreason_go.sh                                      # GO data (~15MB)
```

needs the `mmseqs` cli on PATH (`brew install mmseqs2`).

## build

```bash
.venv/bin/python src/build_train_set.py   # writes the four align*.csv
```

for each benchmark: collect its test seqs, mmseqs-search qfit against them
(sensitive, conservative), drop any qfit seq with a >20% identity hit.

## use

```python
import sys; sys.path.insert(0, "src")
from datasets import QfitDataset, DdgDataset, DmsDataset, GoDataset

align = QfitDataset()                                  # {sequence, uniprot, pdb, chain}

ddg = DdgDataset(dataset="megascale")                  # filter optional
# {sequence, mutation, ddg, pdb, dataset}

dms = DmsDataset(assay="Binding")                      # or dms_id="...", filters optional
# {sequence, mutation, score, assay, dms_id, uniprot}

go = GoDataset()                                       # temporal-holdout test set
# {sequence, protein_id, organism, go_bp, go_mf, go_cc}
```

`QfitDataset(path=...)` to pick a different align set. `DmsDataset` assay is one
of Activity|Binding|Expression|OrganismalFitness|Stability — filter it (or a
`dms_id`), loading all ~2.7M rows is heavy. `GoDataset` is the GO temporal
holdout eval set; go_bp/go_mf/go_cc are GO ids per ontology aspect.
