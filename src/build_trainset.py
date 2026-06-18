# build a leakage-free training set: sample uniref90, drop anything that leaks a functional test
# set (aligns over >COV of a test seq at >MIN_ID identity), and resample until the quota is met.
# checking a sample against the small test sets is cheap, so uniref90 is never scanned whole.
# writes the final set as a fasta keyed by uniprot code.
import random, subprocess, sys, tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as c


def mmseqs(*args):  # run an mmseqs subcommand quietly
    subprocess.run([str(c.mmseqs), *args, "--threads", c.threads, "-v", "1"], check=True)


def test_target(tmp):  # build the combined test-set db once (gpu-padded if GPU=1)
    fa = tmp / "tests.fasta"
    fa.write_text("".join(p.read_text() for p in sorted(c.formatted.glob("*.fasta"))))
    db = tmp / "tests_db"
    mmseqs("createdb", str(fa), str(db), "--dbtype", "1")
    if c.gpu:
        mmseqs("makepaddedseqdb", str(db), str(tmp / "tests_pad"))
        return tmp / "tests_pad"
    return db


def leakers(pool, target, tmp):  # ids in pool [(id,seq)] that cover >COV of a test seq at >MIN_ID
    with open(tmp / "pool.fasta", "w") as o:
        for a, s in pool: o.write(f">{a}\n{s}\n")
    qdb, res, m8 = tmp / "pool_db", tmp / "pool_res", tmp / "pool.m8"
    mmseqs("createdb", str(tmp / "pool.fasta"), str(qdb), "--dbtype", "1")
    mmseqs("search", str(qdb), str(target), str(res), str(tmp / "s"),
           "-s", c.sens, "-e", c.evalue, "-c", c.cov, "--cov-mode", c.cov_mode,
           "--max-seqs", "50", *(("--gpu", "1") if c.gpu else ()))
    mmseqs("convertalis", str(qdb), str(target), str(res), str(m8), "--format-output", "query,target,fident")
    return {ln.split("\t")[0] for ln in open(m8) if float(ln.split("\t")[2]) > c.min_id}


def reservoir(n, exclude, rng):  # uniform-random n in-range uniref seqs not already used; one pass over the fasta
    res, seen = [], 0
    for a, s in c.iter_fasta(c.uniref_fa):
        if a in exclude or not (c.min_len <= len(s) <= c.max_len): continue
        seen += 1
        if len(res) < n:
            res.append((a, s))
        else:
            j = rng.randrange(seen)
            if j < n: res[j] = (a, s)
        if seen % 20_000_000 == 0: print(f"    ...scanned {seen:,}", flush=True)
    return res


def code(uid):  # uniref90 cluster id -> representative uniprot accession
    return uid.split("_", 1)[1] if uid.startswith("UniRef90_") else uid


if __name__ == "__main__":
    rng = random.Random(c.seed)
    c.work.mkdir(parents=True, exist_ok=True)
    c.trainset.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=c.work) as tmp:
        tmp = Path(tmp)
        target = test_target(tmp)
        kept, used, rnd = {}, set(), 0          # code->seq, sampled uniref ids, round counter
        while len(kept) < c.quota:
            rnd += 1
            n = int((c.quota - len(kept)) * c.oversample)
            print(f"round {rnd}: sampling {n:,} candidates (have {len(kept):,}/{c.quota:,})", flush=True)
            pool = reservoir(n, used, rng)
            if not pool:
                print("uniref90 exhausted", flush=True); break
            used.update(a for a, _ in pool)
            bad = leakers(pool, target, tmp)
            print(f"  dropped {len(bad):,} leakers of {len(pool):,}", flush=True)
            for a, s in pool:
                if a not in bad and len(kept) < c.quota: kept[code(a)] = s
        with open(c.trainset, "w") as o:
            for cd, s in kept.items(): o.write(f">{cd}\n{s}\n")
    print(f"wrote {c.trainset} — {len(kept):,} leakage-free seqs", flush=True)
