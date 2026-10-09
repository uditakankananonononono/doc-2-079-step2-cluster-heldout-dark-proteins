"""DOC-2-079 step 2: build armA.tsv and armB.tsv plus FASTA files. Needs frozen inputs: dark_raw.tsv (step 1, md5 a9c57ef14d325e9f7bfc1dea650dabca), uniprot_raw.tsv (DOC-2-072, md5 aa086b6e99d5e84c44789b1a85bc71f1), nopfam_raw.tsv."""
import hashlib, io, numpy as np, pandas as pd
for fn, h in [("dark_raw.tsv", "a9c57ef14d325e9f7bfc1dea650dabca"), ("uniprot_raw.tsv", "aa086b6e99d5e84c44789b1a85bc71f1")]:
    assert hashlib.md5(open(fn).read().encode()).hexdigest() == h, fn
bad = "[XBZUO]"
def one_pfam(d, pc):
    d["pf"] = d[pc].fillna("").str.strip(";").str.split(";"); d = d[d.pf.str.len() == 1].copy(); d["family"] = d.pf.str[0]; return d
c = pd.read_csv("uniprot_raw.tsv", sep="\t"); c.columns = ["accession", "ec", "pfam", "length", "sequence"]; c = one_pfam(c, "pfam"); c["group"] = "characterized"
k = pd.read_csv("dark_raw.tsv", sep="\t"); k.columns = ["accession", "name", "pfam", "length", "sequence"]; k = one_pfam(k, "pfam"); k["group"] = "dark"
n = pd.read_csv("nopfam_raw.tsv", sep="\t"); n.columns = ["accession", "name", "length", "sequence"]; n["family"] = ""; n["group"] = "dark_nopfam"
cols = ["accession", "family", "group", "length", "sequence"]
for x in (c, k, n): x.drop(x[x.sequence.str.contains(bad)].index, inplace=True)
c = c[~c.accession.isin(set(k.accession) | set(n.accession))]
# arm A: dedup exact sequences (keep first accession), families with >=5 dark and >=5 characterized, <=20 each per family (seed 801)
a = pd.concat([c[cols], k[cols]]).sort_values("accession").drop_duplicates("sequence")
cnt = a.groupby(["family", "group"]).size().unstack(fill_value=0); ok = sorted(cnt[(cnt["dark"] >= 5) & (cnt["characterized"] >= 5)].index)
rng = np.random.default_rng(801); out = []
for f in ok:
    for g in ("dark", "characterized"):
        s = a[(a.family == f) & (a.group == g)]; out.append(s.iloc[np.sort(rng.choice(len(s), min(20, len(s)), replace=False))])
A = pd.concat(out); A.to_csv("armA.tsv", sep="\t", index=False)
# arm B: 1,500 dark_nopfam + 1,500 characterized (EC, any Pfam), exact-sequence dedup across both, seed 802
b = pd.concat([n[cols], c[cols]]).sort_values("accession").drop_duplicates("sequence"); rng = np.random.default_rng(802); out = []
for g, m in (("dark_nopfam", 1500), ("characterized", 1500)):
    s = b[b.group == g]; out.append(s.iloc[np.sort(rng.choice(len(s), min(m, len(s)), replace=False))])
B = pd.concat(out); B.to_csv("armB.tsv", sep="\t", index=False)
for nm, D in (("armA", A), ("armB", B)):
    with open(f"{nm}.fasta", "w") as fh:
        for r in D.itertuples(): fh.write(f">{r.accession}\n{r.sequence}\n")
print("BUILD_DONE armA families", len(ok), "n", len(A), A.group.value_counts().to_dict(), "| armB n", len(B), B.group.value_counts().to_dict())
