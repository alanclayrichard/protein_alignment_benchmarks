# protein-alignment-benchmarks

build a **leakage-free** subset of uniref90 to align/train a protein sequence model
on, then test whether alignment improves prediction on functional benchmarks
(thermompnn ddG, proteingym dms fitness, bioreason-pro GO function).

leakage rule: a uniref90 seq leaks a benchmark if one of its **test** seqs aligns
over >80% of that test seq's length at >20% identity (`COV`/`MIN_ID`,
`--cov-mode 2`) — i.e. the test protein's content is present in the uniref seq, so
training on it would leak the eval. we scan test seqs only (eval is test-only).

## pipeline

```
data/download_*.sh    download a benchmark's test set            -> data/
format_benchmark.py   a benchmark's test set -> uniform fasta    -> data/formatted/<bench>.fasta
scan_leakage.py       <bench>.fasta vs uniref90 -> leaked ids     -> leakage/<bench>.txt
build_trainset.py     uniref90 minus leaked ids, sampled         -> trainset.fasta
datasets.py           load trainset + benchmarks as torch datasets
```

uniref90 is never copied: the leakage lists are just ids, and the training set is
whatever sampled size you ask for.

## setup

uniref90 is downloaded separately — point `.env` at it (and at the mmseqs binary).
copy `.env.example` to `.env` and edit; `config.py` auto-loads it. then:

```bash
uv sync                                  # create .venv with core deps (add --extra torch for datasets.py)
source .venv/bin/activate                # then plain `python ...` uses the venv
sh data/download_ddg.sh && sh data/download_proteingym.sh && sh data/download_bioreason_go.sh
```

mmseqs2 isn't pip-installable — set `MMSEQS` to its path, or `brew install mmseqs2`.

## run

```bash
python src/format_benchmark.py           # all benchmarks -> data/formatted/*.fasta
python src/scan_leakage.py               # scan each vs uniref90 -> leakage/*.txt   (the slow step)
python src/build_trainset.py             # uniref90 - leaked ids, sampled -> trainset.fasta
```

sampling lives in `.env`: `TRAIN_SAMPLE=0.05` keeps a random 5%, `TRAIN_SIZE=1000000`
caps the count — so the training set is as big as you want without pulling all ~121M.


## use

needs torch (`uv sync --extra torch`), with the venv activated:

```python
import sys; sys.path.insert(0, "src")
from datasets import TrainSet, DdgDataset, DmsDataset, GoDataset

train = TrainSet()                         # {sequence, id} — the leakage-free align set
ddg   = DdgDataset(dataset="megascale")    # {sequence, mutation, ddg, pdb, dataset}
dms   = DmsDataset(assay="Binding")        # {sequence, mutation, score, assay, dms_id, uniprot}
go    = GoDataset()                        # {sequence, protein_id, organism, go_bp, go_mf, go_cc}
```

`DmsDataset` filters by assay (Activity|Binding|Expression|OrganismalFitness|Stability)
or `dms_id` — loading all ~2.7M rows is heavy.

## datasets

| benchmark | source | unique test seqs | median len |
|---|---|---|---|
| ddg | thermompnn megascale + fireprot | 28,227 | 57 aa |
| dms | proteingym DMS_substitutions | 187 | 222 aa |
| go  | bioreason-pro (cafa5 temporal holdout) | 8,528 | 588 aa |

align pool: uniref90, 121,389,642 seqs. 

## adding a benchmark

1. add `data/download_<name>.sh` to fetch its test set.
2. add an entry to `ADAPTERS` in `src/format_benchmark.py` (name -> fn yielding its
   test seqs), then `python src/format_benchmark.py <name>`.
3. `python src/scan_leakage.py` — only the new benchmark is scanned (the rest are
   cached; `FORCE=1` rescans).
4. `python src/build_trainset.py` — rebuilds the training set leakage-free against
   every benchmark, now including the new one.

## running the scan at scale

`scan_leakage.py` is CPU-bound mmseqs over ~121M sequences, so run it on a machine
with plenty of cores + RAM (16+ cores helps a lot). set `THREADS`/`SPLIT_MEM` in
`.env` to fit the box, point `UNIREF_*` at your uniref90, and detach it:

```bash
nohup python src/scan_leakage.py > "$LEAKAGE_DIR/scan.log" 2>&1 &
```

it's resumable — each benchmark's list is written as it finishes, and a re-run
skips lists that already exist.
