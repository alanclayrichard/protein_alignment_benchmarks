# build a leakage-free training set: stream uniref90, drop every leaked id, sample to size.
# safe = uniref90 minus the union of the chosen leakage lists (default: all in <leakage>).
# TRAIN_SAMPLE keeps a random fraction and TRAIN_SIZE caps the count, so the set is as big as
# you want without pulling all ~121M. pass benchmark names to use just those leakage lists.
# writes <trainset> fasta. re-run after adding a benchmark to refresh the set.
import random, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as c


def leaked_ids(names):  # (benchmark names used, set of leaked uniref ids)
    lists = [c.leakage / f"{n}.txt" for n in names] if names else sorted(c.leakage.glob("*.txt"))
    leaked = set()
    for p in lists:
        leaked |= {l.strip() for l in open(p) if l.strip()}
    return [p.stem for p in lists], leaked


if __name__ == "__main__":
    names, leaked = leaked_ids(sys.argv[1:])
    print(f"leakage-free against: {', '.join(names) or '(no lists found!)'} "
          f"— {len(leaked):,} leaked ids excluded", flush=True)
    rng = random.Random(c.seed)
    c.trainset.parent.mkdir(parents=True, exist_ok=True)
    kept = seen = 0
    with open(c.trainset, "w") as o:
        for acc, seq in c.iter_fasta(c.uniref_fa):
            seen += 1
            if acc in leaked: continue
            if c.train_sample < 1.0 and rng.random() >= c.train_sample: continue
            o.write(f">{acc}\n{seq}\n"); kept += 1
            if kept % 1_000_000 == 0: print(f"  ...{kept:,} kept", flush=True)
            if c.train_size and kept >= c.train_size: break
    print(f"wrote {c.trainset} — {kept:,} safe seqs (scanned {seen:,})", flush=True)
