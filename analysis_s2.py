"""DOC-2-079 step 2 frozen analysis (PROTOCOL.md lock-1). Run once. Needs armA.tsv, armA_emb.npy, armA_clu_cluster.tsv, armB.tsv, armB_emb.npy, armB_clu_cluster.tsv, armB_hits.m8."""
import json, numpy as np, pandas as pd
rng = np.random.default_rng(12345); AA = "ACDEFGHIKLMNPQRSTVWY"
def comp(d):
    C = np.array([[s.count(a) / len(s) for a in AA] + [np.log(len(s))] for s in d.sequence]); return (C - C.mean(0)) / (C.std(0) + 1e-9)
def unit(X): return X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-12)
def clusters(d, fn):
    m = pd.read_csv(fn, sep="\t", header=None, names=["rep", "mem"]); mp = dict(zip(m.mem, m.rep)); return np.array([mp[a] for a in d.accession])
def boot(units, fn, nb=2000):
    """units: array of cluster/family ids per protein. Returns list of fn(idx) over resampled units."""
    uu = sorted(set(units)); members = {u: np.where(units == u)[0] for u in uu}; out = []
    for _ in range(nb):
        idx = np.concatenate([members[u] for u in rng.choice(uu, len(uu))]); out.append(fn(idx))
    return out
ci = lambda v: [float(x) for x in np.nanpercentile(v, [2.5, 97.5])]
def decide(gap, bs, ctrl_pass, label_pre=""):
    se = float(np.nanstd(bs)); lo, hi = ci(bs); half = (hi - lo) / 2; mdg = 2.8 * se
    if not ctrl_pass: lab = "INVALID"
    elif half > 0.05: lab = "INCONCLUSIVE-POWER"
    elif gap >= 0.05 and lo > 0: lab = "DARK-IS-DARKER"
    else: lab = "NO-DEFICIT-ABOVE-MDG"
    return dict(gap_char_minus_dark=float(gap), ci=[lo, hi], ci_half_width=float(half), se=se, mdg_80pct_power=float(mdg), label=lab)
res = {}
# ---- Arm A: same Pfam families, leave-cluster-out 1-NN family recovery
A = pd.read_csv("armA.tsv", sep="\t"); EA = np.load("armA_emb.npy"); cl = clusters(A, "armA_clu_cluster.tsv"); fam = A.family.values; dark = (A.group == "dark").values
def lco_hit(X, labels, cl):
    S = unit(X) @ unit(X).T; same = cl[:, None] == cl[None, :]; S[same] = -np.inf; j = S.argmax(1)
    ev = np.array([((labels == labels[i]) & (cl != cl[i])).any() for i in range(len(labels))]); return (labels[j] == labels) & ev, ev
hitE, ev = lco_hit(EA, fam, cl); hitC, _ = lco_hit(comp(A), fam, cl)
perm = rng.permutation(fam); hitP, evP = lco_hit(EA, perm, cl)
nf = len(set(fam)); accf = lambda h, m, e: float(h[m & e].sum() / max((m & e).sum(), 1))
a_ctrl = accf(hitP, dark, evP) < 3 / nf and accf(hitP, ~dark, evP) < 3 / nf
gapA = accf(hitE, ~dark, ev) - accf(hitE, dark, ev)
bsA = boot(fam, lambda i: accf(hitE[i], ~dark[i], ev[i]) - accf(hitE[i], dark[i], ev[i]))
dE_C = accf(hitE, dark, ev) - accf(hitC, dark, ev); bsEC = boot(fam, lambda i: accf(hitE[i], dark[i], ev[i]) - accf(hitC[i], dark[i], ev[i]))
res["armA"] = dict(n=len(A), n_families=nf, n_clusters=int(len(set(cl))), evaluable={"dark": int((dark & ev).sum()), "characterized": int((~dark & ev).sum())},
    acc={"E_dark": accf(hitE, dark, ev), "E_char": accf(hitE, ~dark, ev), "C_dark": accf(hitC, dark, ev), "C_char": accf(hitC, ~dark, ev)},
    G1=dict(perm_dark=accf(hitP, dark, evP), perm_char=accf(hitP, ~dark, evP), threshold=3 / nf, pass_=bool(a_ctrl)),
    G2_reported=dict(E_minus_C_dark=float(dE_C), ci=ci(bsEC)), gap=decide(gapA, bsA, a_ctrl))
# ---- Arm B: no-Pfam dark vs characterized, ground truth = alignment-defined remote homologs (20-40% identity, e<=1e-3)
B = pd.read_csv("armB.tsv", sep="\t"); EB = np.load("armB_emb.npy"); clB = clusters(B, "armB_clu_cluster.tsv"); grp = (B.group == "dark_nopfam").values
idx = {a: i for i, a in enumerate(B.accession)}; H = pd.read_csv("armB_hits.m8", sep="\t", header=None, names=["q", "t", "id", "e"]); H = H[H.q != H.t]
n = len(B); close = np.zeros((n, n), bool); remote = np.zeros((n, n), bool)
for q, t, f, e in H.itertuples(index=False):
    i, j = idx[q], idx[t]
    if f >= 0.40: close[i, j] = True
    elif f >= 0.20 and e <= 1e-3: remote[i, j] = True
remote &= ~close; evB = remote.any(1)
def hitsB(X):
    S = unit(X) @ unit(X).T; S[np.eye(n, dtype=bool)] = -np.inf; S[close] = -np.inf; S[clB[:, None] == clB[None, :]] = -np.inf; j = S.argmax(1); return remote[np.arange(n), j] & evB
hB = hitsB(EB); hBC = hitsB(comp(B))
rc = rng.integers(0, n, n); cand = ~(close | (clB[:, None] == clB[None, :]) | np.eye(n, dtype=bool))
ctrl = np.array([remote[i, rng.choice(np.where(cand[i])[0])] if cand[i].any() else False for i in range(n)]) & evB
accb = lambda h, m, e: float(h[m & e].sum() / max((m & e).sum(), 1))
b_ctrl = accb(ctrl, grp, evB) < 0.10 and accb(ctrl, ~grp, evB) < 0.10
gapB = accb(hB, ~grp, evB) - accb(hB, grp, evB); bsB = boot(clB, lambda i: accb(hB[i], ~grp[i], evB[i]) - accb(hB[i], grp[i], evB[i]))
res["armB"] = dict(n=n, n_clusters=int(len(set(clB))), evaluable={"dark_nopfam": int((grp & evB).sum()), "characterized": int((~grp & evB).sum())},
    acc={"E_dark_nopfam": accb(hB, grp, evB), "E_char": accb(hB, ~grp, evB), "C_dark_nopfam": accb(hBC, grp, evB), "C_char": accb(hBC, ~grp, evB)},
    G1=dict(random_neighbour_dark=accb(ctrl, grp, evB), random_neighbour_char=accb(ctrl, ~grp, evB), threshold=0.10, pass_=bool(b_ctrl)), gap=decide(gapB, bsB, b_ctrl))
res["LABELS"] = {"armA": res["armA"]["gap"]["label"], "armB": res["armB"]["gap"]["label"]}
print("RESULT_JSON", json.dumps(res, default=float)); open("results.json", "w").write(json.dumps(res, default=float, indent=1))
