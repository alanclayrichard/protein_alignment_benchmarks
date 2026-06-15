# paths + knobs, all env-overridable (see .env.example). nothing personal hardcoded.
import gzip, os
from pathlib import Path

repo = Path(__file__).resolve().parent.parent
data = repo / "data"

# auto-load repo/.env (KEY=VALUE, $REFS + inline #comments ok) so scripts just run
_env = repo / ".env"
if _env.exists():
    for line in _env.read_text().splitlines():
        line = line.split(" #", 1)[0].strip()
        if "=" in line and not line.startswith("#"):
            k, v = line.split("=", 1)
            os.environ.setdefault(k.strip(), os.path.expandvars(v.strip().strip("\"'")))

# uniref90 (downloaded separately) + mmseqs, plus the repo's inputs/outputs
uniref_dir = Path(os.environ.get("UNIREF_DIR", repo / "uniref"))
mmseqs     = os.environ.get("MMSEQS", "mmseqs")
uniref_db  = Path(os.environ.get("UNIREF_DB", uniref_dir / "uniref90_db"))
uniref_fa  = Path(os.environ.get("UNIREF_FASTA", uniref_dir / "uniref90.fasta.gz"))
formatted  = Path(os.environ.get("FORMATTED_DIR", data / "formatted"))   # formatted test fastas
leakage    = Path(os.environ.get("LEAKAGE_DIR", repo / "leakage"))       # <bench>.txt leaked-id lists
trainset   = Path(os.environ.get("TRAINSET", repo / "trainset.fasta"))   # generated leakage-free set
work       = Path(os.environ.get("WORK", repo / "work"))                 # mmseqs scratch (ephemeral)

# mmseqs + leakage knobs
threads   = os.environ.get("THREADS", "4")
maxseqs   = os.environ.get("MAXSEQS", "2000")     # uniref hits kept per test seq
split_mem = os.environ.get("SPLIT_MEM", "12G")    # mmseqs ram budget
sens      = os.environ.get("SENS", "7.5")         # mmseqs sensitivity
evalue    = os.environ.get("EVALUE", "1e-3")      # mmseqs e-value cutoff
cov       = os.environ.get("COV", "0.8")          # min alignment coverage
cov_mode  = os.environ.get("COV_MODE", "2")       # 2 = fraction of the test seq covered
min_id    = float(os.environ.get("MIN_ID", "0.2"))  # >this id to a test seq = leakage

# training-set sampling (so you don't pull all ~121M)
train_sample = float(os.environ.get("TRAIN_SAMPLE", "1.0"))  # fraction of safe seqs to keep
train_size   = int(os.environ.get("TRAIN_SIZE", "0"))        # cap on kept seqs, 0 = no cap
seed         = int(os.environ.get("SEED", "0"))


def iter_fasta(path):  # (id, sequence) pairs from a (optionally gz) fasta
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        acc, seq = None, []
        for line in f:
            if line.startswith(">"):
                if acc is not None: yield acc, "".join(seq)
                acc, seq = line[1:].split(None, 1)[0], []
            else: seq.append(line.strip())
        if acc is not None: yield acc, "".join(seq)
