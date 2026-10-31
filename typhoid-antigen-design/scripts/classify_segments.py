import collections, numpy as np, pandas as pd
IDENT_SHARED, CONS_MIN = 0.8, 0.8
GRP = {"typhi":"typhi","paratyphiA":"paratyphiA","enteritidis":"nts","typhimurium":"nts","enterica":"nts",
       "ecoli":"conf","shigella":"conf","klebsiella":"conf","pseudomonas":"conf"}
ALN = {"hlyE":"aligned_hlye.aln","yncE":"yncE_aligned.aln","cdtB":"cdtB_aligned.aln"}
full = {r.split("\n")[0].replace("_full",""): "".join(r.split("\n")[1:]) for r in open("full_length.faa").read().split(">")[1:]}
def read_aln(f):
    d, cur = {}, None
    for l in open(f):
        l = l.rstrip()
        if l.startswith(">"): cur = l[1:]; d[cur] = []
        else: d[cur].append(l)
    return {k: "".join(v) for k, v in d.items()}
mm = lambda a, b: sum(x != y for x, y in zip(a, b))
lo = lambda L, ref: min((mm(x, ref) for x in L), default=None)
cache = {}
for p, f in ALN.items():
    aln = read_aln(f); ref = full[p]
    refrow = next(s for n, s in aln.items() if n.startswith("typhi") and s.replace("-", "") == ref)
    cache[p] = (aln, refrow, [i for i, c in enumerate(refrow) if c != "-"])
df = pd.read_csv("epitope_ranking_v3.csv")
rows = []
for r in df.itertuples():
    aln, refrow, cols = cache[r.protein]
    sc = cols[r.start - 1:r.end]; seg = lambda s: "".join(s[i] for i in sc)
    refseg = seg(refrow); by = collections.defaultdict(list)
    for n, s in aln.items():
        g = GRP.get(n.split("_")[0])
        if g: by[g].append(seg(s))
    pa = by["paratyphiA"]; cons = cn = ident = None; cand = False
    if pa:
        cons, cn = collections.Counter(pa).most_common(1)[0]
        ident = float(np.mean([x == refseg for x in pa]))
        cls = "pan-typhoidal" if ident >= IDENT_SHARED else ("Typhi-selective" if ident == 0 else "mixed")
        cand = cls == "Typhi-selective" and cn / len(pa) >= CONS_MIN
    else: cls = "no Paratyphi A genomes"
    rows.append(dict(rank_no_prior=r.rank_no_prior, protein=r.protein, start=r.start, end=r.end, length=r.length,
        typhi_peptide=refseg, n_parA=len(pa), parA_identical_frac=None if ident is None else round(ident, 2),
        parA_consensus_peptide=cons, parA_consensus_frac=None if not pa else round(cn / len(pa), 2),
        mm_typhi_vs_parA=None if cons is None else mm(refseg, cons), typhi_class=cls,
        parA_allele_candidate=cand, typhi_vs_nts_min_mm=lo(by["nts"], refseg), typhi_vs_conf_min_mm=lo(by["conf"], refseg),
        parA_allele_vs_nts_min_mm=lo(by["nts"], cons) if cand else None,
        parA_allele_vs_conf_min_mm=lo(by["conf"], cons) if cand else None))
out = pd.DataFrame(rows).sort_values("rank_no_prior"); out.to_csv("epitope_classes.csv", index=False)
print(out[["rank_no_prior","protein","start","end","typhi_class","parA_identical_frac","mm_typhi_vs_parA","parA_allele_candidate","parA_consensus_peptide"]].to_string(index=False))
print("\nclass counts:", out.typhi_class.value_counts().to_dict(), "| Paratyphi A allele candidates:", int(out.parA_allele_candidate.sum()))
