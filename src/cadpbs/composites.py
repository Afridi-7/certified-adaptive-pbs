"""
Composite specimen generator (exact ground truth) for Part I / Part III.

A SPECIMEN is a finite population of N fields. A FIELD is a bag of m cells drawn from KU-Optofil.
  - rarity k   : number of CRITICAL fields (fields containing >= 1 blast) out of N
  - dilution m : cells per field. m=1 means the field IS a cell; larger m dilutes a blast inside its field.
Cells are drawn WITHOUT replacement within a specimen, from ONE cross-fitting fold, so the embeddings and
classifier scores of every cell come from a model that never saw that cell's patient.

Cell pool (pre-registration decision, primary specification): "pool" = background cells drawn from all
blast-free patients of the fold. Sensitivity: "host" = background from a single randomly chosen patient
(with replacement if that patient has too few cells).

Field representation (pre-registration decision, primary specification): MEAN of the cell embeddings.
Sensitivity: "max" (element-wise) and "single" (the cell embedding nearest the mean).
This choice matters: mean-pooling is what dilutes a blast, see docs/separability_v1.md.

Everything needed to audit an episode is returned: the cell indices of every field, the critical mask,
the true blast count per field and the specimen-level blast proportion.
"""
import numpy as np


class CellBank:
    def __init__(self, emb, df, l2_normalise=True):
        self.X = emb / np.linalg.norm(emb, axis=1, keepdims=True) if l2_normalise else emb
        self.df = df.reset_index(drop=True)
        self.is_blast = self.df.is_blast.values.astype(bool)
        self.fold = self.df.fold.values
        self.pid = self.df.patient_id.values
        self.p = self.df.probe_p_oof.values            # cross-fitted classifier score (the teacher signal)
        blast_pat = set(self.df.loc[self.is_blast, "patient_id"])
        self.blast_idx = {f: np.flatnonzero((self.fold == f) & self.is_blast) for f in np.unique(self.fold)}
        self.bg_idx = {f: np.flatnonzero((self.fold == f) & ~self.is_blast &
                                         ~np.isin(self.pid, list(blast_pat))) for f in np.unique(self.fold)}
        self.bg_by_patient = {}
        for f in np.unique(self.fold):
            pids = np.unique(self.pid[self.bg_idx[f]])
            self.bg_by_patient[f] = {p: np.flatnonzero((self.pid == p) & ~self.is_blast) for p in pids}


def make_specimen(bank, fold, N, k, m, rng, pool="pool", pooling="mean", blasts_per_critical=1):
    """Returns dict with field embeddings F (N x d), critical mask, cell indices, blast counts."""
    need_bg = N * m            # draw a full background, then overwrite one cell per critical field
    if pool == "pool":
        src = bank.bg_idx[fold]
        bg = rng.choice(src, need_bg, replace=need_bg > len(src))
    elif pool == "host":
        pats = list(bank.bg_by_patient[fold]); host = pats[rng.integers(len(pats))]
        src = bank.bg_by_patient[fold][host]
        bg = rng.choice(src, need_bg, replace=need_bg > len(src))
    else:
        raise ValueError(pool)
    bl = rng.choice(bank.blast_idx[fold], k * blasts_per_critical,
                    replace=k * blasts_per_critical > len(bank.blast_idx[fold]))
    cells = np.asarray(bg, int).reshape(N, m)
    crit = np.zeros(N, bool)
    crit_fields = rng.choice(N, k, replace=False)
    crit[crit_fields] = True
    for j, fi in enumerate(crit_fields):                # overwrite one cell of each critical field
        for b in range(blasts_per_critical):
            cells[fi, b] = bl[j * blasts_per_critical + b]
    Xc = bank.X[cells]                                   # (N, m, d)
    if pooling == "mean":
        F = Xc.mean(1)
    elif pooling == "max":
        F = Xc.max(1)
    elif pooling == "single":
        mu = Xc.mean(1, keepdims=True)
        F = Xc[np.arange(N), np.argmax((Xc * mu).sum(-1), axis=1)]
    else:
        raise ValueError(pooling)
    F = F / np.linalg.norm(F, axis=1, keepdims=True)
    nb = bank.is_blast[cells].sum(1)
    teacher = bank.p[cells].max(1)                       # field relevance from the cross-fitted classifier
    return dict(F=F, crit=crit, cells=cells, n_blast=nb, teacher=teacher,
                blast_fraction=float(nb.sum() / (N * m)), fold=fold, N=N, k=k, m=m,
                pool=pool, pooling=pooling)
