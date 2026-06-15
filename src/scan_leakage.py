# scan formatted test sets against uniref90 and write the leakage lists (runs on the pod).
# for each <formatted>/<name>.fasta: mmseqs-search vs uniref90, write the leaked uniref ids
# (a hit covering >COV of a test seq at >MIN_ID identity = leakage) to <leakage>/<name>.txt.
# auto-discovers every formatted test set; re-run skips ones already scanned (FORCE=1 rescans),
# so adding a benchmark only scans the new one. pass names to scan just those.
import os, subprocess, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as c

FORCE = os.environ.get("FORCE", "")


def mmseqs(*args):
    subprocess.run([str(c.mmseqs), *args, "--threads", c.threads, "-v", "2"], check=True)


def scan(fasta):  # write the leaked uniref ids for one formatted test set, return the count
    name = fasta.stem
    out = c.leakage / f"{name}.txt"
    if out.exists() and not FORCE:
        return sum(1 for _ in open(out))
    db, res, m8 = (c.work / f"{name}{x}" for x in ("_db", "_res", ".m8"))
    mmseqs("createdb", str(fasta), str(db), "--dbtype", "1")
    mmseqs("search", str(db), str(c.uniref_db), str(res), str(c.work / f"{name}_tmp"),
        "-s", c.sens, "-e", c.evalue, "-c", c.cov, "--cov-mode", c.cov_mode,
        "--max-seqs", c.maxseqs, "--split-memory-limit", c.split_mem)
    mmseqs("convertalis", str(db), str(c.uniref_db), str(res), str(m8),
        "--format-output", "query,target,fident")
    # collect leaked uniref ids and count per-query hits (to spot --max-seqs saturation)
    leaked, per_q = set(), {}
    for line in open(m8):
        q, t, fid = line.split("\t")[:3]
        per_q[q] = per_q.get(q, 0) + 1
        if float(fid) > c.min_id: leaked.add(t)
    capped = sum(n >= int(c.maxseqs) for n in per_q.values())
    if capped: print(f"  {name}: {capped} queries hit the --max-seqs cap "
                     f"({c.maxseqs}) — bump MAXSEQS to be sure", flush=True)
    out.write_text("\n".join(sorted(leaked)))
    return len(leaked)


if __name__ == "__main__":
    c.leakage.mkdir(parents=True, exist_ok=True)
    c.work.mkdir(parents=True, exist_ok=True)
    fastas = [c.formatted / f"{n}.fasta" for n in sys.argv[1:]] or sorted(c.formatted.glob("*.fasta"))
    if not fastas:
        sys.exit("no formatted test sets in %s — run format_benchmark.py first" % c.formatted)
    for fa in fastas:
        n = scan(fa)
        print(f"{fa.stem}: {n} leaked uniref ids -> {c.leakage}/{fa.stem}.txt", flush=True)
