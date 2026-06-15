# torch datasets: the leakage-free training set + the ddG / dms / GO eval benchmarks.
import ast, csv, sys
from pathlib import Path
import pandas as pd
from torch.utils.data import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as c

csv.field_size_limit(sys.maxsize)
repo = c.repo


class TrainSet(Dataset):
    # the generated leakage-free training set (make it with build_trainset.py).
    # item: sequence, id. it's already sampled to size, so we load it in memory.
    def __init__(self, path=c.trainset):
        self.items = [{"sequence": s, "id": a} for a, s in c.iter_fasta(path)]

    def __len__(self): return len(self.items)

    def __getitem__(self, i): return self.items[i]


class DdgDataset(Dataset):
    # thermompnn ddG test set. item: sequence, mutation, ddg, pdb, dataset.
    # dataset is megascale|fireprot (pass to filter); sign follows thermompnn.
    files = {
        "megascale": "testing/mega_test.csv",
        "fireprot":  "testing/fireprot_test.csv",
    }

    def __init__(self, dataset=None, data_dir=repo / "data" / "ThermoMPNN" / "data_all"):
        self.items = []
        for ds, rel in self.files.items():
            if dataset and ds != dataset: continue
            self.items += list(self._load(Path(data_dir) / rel, ds))

    def _load(self, path, ds):
        if not path.exists(): return
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                it = self._megascale(row) if ds == "megascale" else self._fireprot(row)
                if it: it["dataset"] = ds; yield it

    @staticmethod
    def _megascale(r):
        mut, ddg = r.get("mut_type", "").strip(), r.get("ddG_ML", "").strip()
        if not mut or mut == "wt" or ddg in ("", "-"): return None
        if any(x in mut for x in ("ins", "del", ":")): return None  # single subs only
        return {"sequence": r["wt_seq"].strip().upper(), "mutation": mut,
                "ddg": -float(ddg), "pdb": r.get("WT_name", "").replace(".pdb", "")}

    @staticmethod
    def _fireprot(r):
        if not (r.get("ddG") or "").strip(): return None
        return {"sequence": r["pdb_sequence"].strip().upper(),
                "mutation": f"{r['wild_type']}{r['position']}{r['mutation']}",
                "ddg": float(r["ddG"]), "pdb": r.get("pdb_id_corrected", "")}

    def __len__(self): return len(self.items)

    def __getitem__(self, i): return self.items[i]


class GoDataset(Dataset):
    # bioreason-pro GO function test set (cafa5, temporal holdout).
    # item: sequence, protein_id, organism, go_bp, go_mf, go_cc (lists of GO ids).
    def __init__(self, data_dir=repo / "data" / "bioreason_go"):
        df = pd.read_parquet(Path(data_dir) / "test" / "data" / "test-00000-of-00001.parquet")
        self.items = [{"sequence": str(r.sequence).strip().upper(),
                       "protein_id": r.protein_id, "organism": r.organism,
                       "go_bp": _terms(r.go_bp), "go_mf": _terms(r.go_mf),
                       "go_cc": _terms(r.go_cc)}
                      for r in df.itertuples()]

    def __len__(self): return len(self.items)

    def __getitem__(self, i): return self.items[i]


def _terms(v):  # go column: array or list-repr string, nan if absent
    if isinstance(v, str): return ast.literal_eval(v) if v.startswith("[") else []
    return list(v) if hasattr(v, "__len__") else []


class DmsDataset(Dataset):
    # proteingym dms data. item: sequence, mutation, score, assay, dms_id, uniprot.
    # pass assay (Activity|Binding|Expression|OrganismalFitness|Stability) or dms_id to filter.
    def __init__(self, assay=None, dms_id=None, data_dir=repo / "data" / "proteingym"):
        data_dir = Path(data_dir)
        self.items = []
        with open(data_dir / "DMS_substitutions.csv", newline="") as f:
            meta = {r["DMS_id"]: (r["target_seq"].strip().upper(),
                                  r["coarse_selection_type"], r["UniProt_ID"])
                    for r in csv.DictReader(f)}
        for did, (seq, atype, uni) in meta.items():
            if (assay and atype != assay) or (dms_id and did != dms_id): continue
            path = data_dir / "DMS_ProteinGym_substitutions" / f"{did}.csv"
            if not path.exists(): continue
            with open(path, newline="") as f:
                for r in csv.DictReader(f):
                    if not (r.get("DMS_score") or "").strip(): continue
                    self.items.append({"sequence": seq, "mutation": r["mutant"],
                        "score": float(r["DMS_score"]), "assay": atype,
                        "dms_id": did, "uniprot": uni})

    def __len__(self): return len(self.items)

    def __getitem__(self, i): return self.items[i]
