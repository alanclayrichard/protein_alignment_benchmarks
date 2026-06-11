"""build leakage-free qfit align sets: drop qfit seqs >20% identical to a benchmark seq.

writes align_ddg_safe.csv (vs ddG), align_dms_safe.csv (vs proteingym dms),
align_go_safe.csv (vs bioreason-pro GO), and align.csv (safe from all three).
"""
import csv, subprocess, sys
from pathlib import Path
import pandas as pd

csv.field_size_limit(sys.maxsize)
repo = Path(__file__).resolve().parent.parent
ddg_dir = repo / "data" / "ThermoMPNN" / "data_all"
dms_dir = repo / "data" / "proteingym"
go_dir = repo / "data" / "bioreason_go"
work = repo / "work"; work.mkdir(exist_ok=True)


def col(path, *names):  # yield non-empty values from the given columns
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            for n in names:
                v = (row.get(n) or "").strip().upper()
                if v: yield v


def fasta(path, items):
    with open(path, "w") as o:
        for name, seq in items: o.write(f">{name}\n{seq}\n")


# unique qfit seqs (keep first row per seq)
qfit = {}
with open(repo / "data" / "qfit_sequences.csv", newline="") as f:
    for row in csv.DictReader(f):
        qfit.setdefault(row["sequence"].strip().upper(), row)
ids = {s: f"q{i}" for i, s in enumerate(qfit)}
fasta(work / "qfit.fasta", ((ids[s], s) for s in qfit))


def unsafe_against(ref_seqs, name):
    """qfit ids with any >20% identity local hit to ref_seqs (sensitive, conservative)."""
    fasta(work / f"{name}.fasta", ((f"r{i}", s) for i, s in enumerate(ref_seqs)))
    hits = work / f"hits_{name}.m8"
    subprocess.run(["mmseqs", "easy-search", str(work / "qfit.fasta"), str(work / f"{name}.fasta"),
        str(hits), str(work / "tmp"), "-s", "7.5", "-e", "1", "-c", "0",
        "--max-seqs", "4000", "--threads", "8", "-v", "1",
        "--format-output", "query,fident"], check=True, stdout=subprocess.DEVNULL)
    return {l.split("\t")[0] for l in open(hits) if float(l.split("\t")[1]) > 0.2}


# refs are TEST seqs only — we only eval on test, so leakage into train is irrelevant
# (and using train over-prunes qfit, esp. the 133k-seq GO train set).

# ddG ref = thermompnn megascale + fireprot test splits
ddg = set(col(ddg_dir / "testing/mega_test.csv", "aa_seq", "wt_seq"))
ddg |= set(col(ddg_dir / "testing/fireprot_test.csv", "pdb_sequence"))

# dms ref = proteingym target seqs (eval-only benchmark, no train split)
dms = set(col(dms_dir / "DMS_substitutions.csv", "target_seq"))

# go ref = bioreason-pro temporal-holdout test seqs
go = set(s.strip().upper() for s in
         pd.read_parquet(go_dir / "test/data/test-00000-of-00001.parquet", columns=["sequence"])["sequence"] if s)

unsafe_ddg = unsafe_against(ddg, "ddg")
unsafe_dms = unsafe_against(dms, "dms")
unsafe_go = unsafe_against(go, "go")

cols = ["pdb_code", "chain_id", "uniprot_id", "sequence"]


def write(path, unsafe):
    kept = [r for s, r in qfit.items() if ids[s] not in unsafe]
    with open(path, "w", newline="") as o:
        w = csv.DictWriter(o, fieldnames=cols); w.writeheader()
        w.writerows({c: r.get(c, "") for c in cols} for r in kept)
    print(f"{path.name}: kept {len(kept)}/{len(qfit)}")


write(repo / "align_ddg_safe.csv", unsafe_ddg)
write(repo / "align_dms_safe.csv", unsafe_dms)
write(repo / "align_go_safe.csv", unsafe_go)
write(repo / "align.csv", unsafe_ddg | unsafe_dms | unsafe_go)  # safe from all three
