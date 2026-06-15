# format functional-benchmark TEST sets into uniform fastas for leakage scanning.
# each benchmark = one ADAPTERS entry: name -> fn yielding its test seqs (from anywhere).
# writes deduped, uppercased seqs to <formatted>/<name>.fasta.
# add a benchmark: add an entry (+ a data/download_<name>.sh), then
#   python src/format_benchmark.py <name>
import csv, sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as c

csv.field_size_limit(sys.maxsize)


def _col(path, *names):  # non-empty upper values from the given csv columns
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            for n in names:
                v = (row.get(n) or "").strip().upper()
                if v: yield v


def ddg():  # thermompnn megascale + fireprot test splits
    d = c.data / "ThermoMPNN" / "data_all" / "testing"
    yield from _col(d / "mega_test.csv", "aa_seq", "wt_seq")
    yield from _col(d / "fireprot_test.csv", "pdb_sequence")


def dms():  # proteingym target seqs (eval-only benchmark)
    yield from _col(c.data / "proteingym" / "DMS_substitutions.csv", "target_seq")


def go():  # bioreason-pro temporal-holdout test seqs
    p = c.data / "bioreason_go" / "test" / "data" / "test-00000-of-00001.parquet"
    for s in pd.read_parquet(p, columns=["sequence"])["sequence"]:
        if s: yield str(s).strip().upper()


ADAPTERS = {"ddg": ddg, "dms": dms, "go": go}  # name -> test-seq iterator


def write(name):  # dedup + write <formatted>/<name>.fasta, return (path, count)
    c.formatted.mkdir(parents=True, exist_ok=True)
    out = c.formatted / f"{name}.fasta"
    n = 0
    with open(out, "w") as o:
        for i, s in enumerate(dict.fromkeys(ADAPTERS[name]())):
            o.write(f">{name}_{i}\n{s}\n"); n += 1
    return out, n


if __name__ == "__main__":
    for name in sys.argv[1:] or list(ADAPTERS):
        out, n = write(name)
        print(f"{name}: {n} test seqs -> {out}")
