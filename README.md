# protein-alignment-benchmarks

build a **leakage-free** subset of uniref90 to align/train a protein sequence model on, then test whether alignment improves prediction on functional benchmarks (thermompnn ddG, proteingym dms fitness, bioreason-pro GO function, allobench + passerrank allosteric sites, human PPI, and whatever else you want).

leakage rule: a sampled uniref90 sequence is dropped if it aligns over >80% of a **test** seq's length at >20% identity (tunable with `COV`/`MIN_ID`) — i.e. the test protein's content is present in the train seq, so training on it would leak the eval.

## leakage coverage modes (`COV_MODE`)

mmseqs `--cov-mode` decides which sequence the `COV` (0.8) threshold applies to, which matters for **multidomain** proteins:

- `COV_MODE=1` (default) — cover the **test** seq. drops a train seq if it covers >80% of a test protein *regardless of the train protein's length*, so it catches a long multidomain uniref protein that contains a whole test protein as one domain. the leakage-faithful choice.
- `COV_MODE=0` — cover **both** seqs (>80% of each). only removes near-duplicates of similar length; it would **keep** the multidomain protein above (the test domain is <80% of the long train seq), letting a test protein leak in as a sub-domain.
- `COV_MODE=2` — cover the **train** seq. catches the reverse (a short train protein that is a sub-domain of a longer test protein).

for example: a 600 aa uniref protein whose 200 aa domain B equals a 200 aa test protein — the alignment covers 100% of the test but only 33% of the train seq, so `COV_MODE=1` flags it (removed) while `COV_MODE=0` keeps it.

## usage pipeline
to create your own leakage-free training set:

```
data/download_*.sh    download a benchmark's test set            -> data/<bench>/
format_benchmark.py   each test set -> uniform fasta             -> data/formatted/<bench>.fasta
build_trainset.py     sample uniref90, drop leakers, repeat      -> trainset.fasta
datasets.py           load trainset + benchmarks as torch datasets
```

instead of scanning all ~121M uniref90 seqs, `build_trainset.py` optionally **samples** candidates, checks just that sample against the (small) test sets with mmseqs, keeps the clean ones, and iteratively resamples until the quota is met. checking a sample against the test sets is cheap, so it's fast — and exact. uniref90 is never copied; only the sampled training set is written. The sampling strategy is currently uniform sampling over the whole uniref90 to ensure coverage of protein lengths and families but could be optionally tuned to stratify on a certain feature (e.g. sequence composition, family, fold, organism, etc.).

## setup

uniref90 is downloaded separately (https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/) — point `.env` at it (and at the mmseqs binary). copy `.env.example` to `.env` and edit; `config.py` auto-loads it. then:

```bash
uv sync                                  # create .venv with core deps (add --extra torch for datasets.py)
source .venv/bin/activate                # then plain `python ...` uses the venv
for s in data/download_*.sh; do sh "$s"; done   # fetch every benchmark's test set
```

set `MMSEQS` to its path (or `brew install mmseqs2`). for the gpu check, point `MMSEQS` at a gpu build and set `GPU=1`.

## run

```bash
python src/format_benchmark.py           # all benchmarks -> data/formatted/*.fasta
python src/build_trainset.py             # sample uniref90, drop leakers, repeat -> trainset.fasta
```

`.env` controls it: `QUOTA` is how many leakage-free seqs to collect for the training set, `OVERSAMPLE` how much extra to draw each round to cover dropped leakers, `MIN_LEN`/`MAX_LEN` bound the sampled protein length (drops fragments + giant non-physiological seqs), `GPU=1` runs the check on a cuda gpu.

## use

needs torch (`uv sync --extra torch`), with the venv activated:

```python
import sys; sys.path.insert(0, "src")
from datasets import (TrainSet, DdgDataset, DmsDataset, GoDataset,
                      AllobenchDataset, PpiDataset, PasserrankDataset,
                      PdbbindDataset, BindingdbDataset)

train = TrainSet()                         # {sequence, id} — the leakage-free align set (id = uniprot code)
ddg   = DdgDataset(dataset="megascale")    # {sequence, mutation, ddg, pdb, dataset}
dms   = DmsDataset(assay="Binding")        # {sequence, mutation, score, assay, dms_id, uniprot}
go    = GoDataset()                         # {sequence, protein_id, organism, go_bp, go_mf, go_cc}
allo  = AllobenchDataset()                  # {sequence, target_id, gene, organism, uniprot, allosteric_site, active_site}
ppi   = PpiDataset(level=0)                 # {sequence_a, sequence_b, id_a, id_b, label, level}
pr    = PasserrankDataset()                 # {sequence, uniprot, gene, organism, pdb, allosteric_site}
pdb   = PdbbindDataset()                    # {sequence, pdb, value, smiles, split} — LP-PDBBind affinity, test split
bdb   = BindingdbDataset()                  # {sequence, uniprot, target_name, organism, smiles, measure, value, relation}
```

## datasets

| benchmark | source | unique test seqs |
|---|---|---|
| ddg | thermompnn megascale + fireprot | 28,227 |
| dms | proteingym DMS_substitutions | 187 |
| go  | bioreason-pro (cafa5 temporal holdout) | 8,528 |
| allobench | allobench allosteric/active sites | 425 |
| ppi | human ppi gold standard (figshare) | 11,019 |
| passerrank | passerrank allosteric set (ASD) | 333 |
| pdbbind | LP-PDBBind protein–ligand affinity (test split) | 2,644 |
| bindingdb | bindingdb articles, protein–ligand affinity (single-chain targets) | 2,157 |

sample pool: uniref90, ~121M sequences.

## adding a benchmark

1. add `data/download_<name>.sh` to fetch its test set.
2. add an adapter to `ADAPTERS` in `src/format_benchmark.py` (name -> fn yielding its test seqs),
   then `python src/format_benchmark.py <name>`.
3. `python src/build_trainset.py` — the sample is checked against every `data/formatted/*.fasta`,
   so the new benchmark is protected automatically.


## TODO:
- add other sampling strategies for train set  
- add more benchmarks
- add structual leakge checking with TM-align scoring when structures are available
- add script to get uniprot training entries with qfit reanalyzed structures 
