"""Rank consensus epitope segments using the genome alignments (prevalence-weighted).
Criteria (weights): spec_nts .30, spec_other .20, accessibility .15, conservation .25, prior_evidence .10 (optional).
Weights are renormalised over the criteria available. usage: python3 -I rank_epitopes_v2.py"""
import os, collections, numpy as np, pandas as pd
rng=np.random.default_rng(42)
W={"spec_nts":.30,"spec_other":.20,"accessibility":.15,"conservation":.25,"prior_evidence":.10}
ALN={"hlyE":"aligned_hlye.aln","yncE":"yncE_aligned.aln","cdtB":"cdtB_aligned.aln"}
GRP={"typhi":"typhi","paratyphiA":"paratyphiA","enteritidis":"nts","typhimurium":"nts","enterica":"nts","ecoli":"conf","shigella":"conf","klebsiella":"conf"}
full={r.split("\n")[0].replace("_full",""):"".join(r.split("\n")[1:]) for r in open("full_length.faa").read().split(">")[1:]}
gl=pd.read_csv("gene_level_summary.csv").set_index("gene")
def read_aln(f):
    d={}; cur=None
    for l in open(f):
        l=l.rstrip()
        if l.startswith(">"): cur=l[1:]; d[cur]=[]          # duplicate headers carry identical sequences -> collapse
        else: d[cur].append(l)
    return {k:"".join(v) for k,v in d.items()}
share=lambda mm: max(0.0,1-mm/3)                          # same 0-3 mismatch scale as the user's S_conf
rows=[]
for r in pd.read_csv("consensus_epitopes.csv").itertuples():
    if r.length<7: continue
    aln=read_aln(ALN[r.protein]); ref=full[r.protein]
    refrow=next(s for n,s in aln.items() if n.startswith("typhi") and s.replace("-","")==ref)
    cols=[i for i,c in enumerate(refrow) if c!="-"]       # alignment column of each reference position (insertions skipped)
    sc=cols[r.start-1:r.end]; seg="".join(refrow[i] for i in sc); assert seg==ref[r.start-1:r.end]
    mm=collections.defaultdict(list)
    for n,s in aln.items(): mm[GRP[n.split("_")[0]]].append(sum(s[i]!=refrow[i] for i in sc))
    g=gl.loc[r.protein]
    carry=lambda grp,pct: pct/100*np.mean([share(m) for m in mm[grp]]) if mm[grp] else 0.0
    p_nts=carry("nts",g.nts_pct); p_conf=carry("conf",g.confounder_pct)
    cons=g.typhoidal_pct/100*np.mean([m==0 for m in mm["typhi"]])
    pa=mm["paratyphiA"]
    rows.append(dict(protein=r.protein,start=r.start,end=r.end,length=r.length,peptide=seg,
        conservation=cons,p_share_nts=p_nts,p_share_conf=p_conf,n_nts=len(mm["nts"]),n_conf=len(mm["conf"]),
        paratyphiA_min_mm=min(pa) if pa else None,
        klass="pan-typhoidal" if pa and min(pa)==0 else "Typhi-selective"))
df=pd.DataFrame(rows)
# non-homologous chance matches (human, Plasmodium, Pseudomonas) from the exact-match screen
cs=pd.read_csv("epitope_crossreactivity_summary.csv")
cs["pen_nonhomolog"]=np.clip((cs[["longest_human","longest_plasmodium","longest_other_bacteria"]].max(axis=1)-6)/4,0,1)
df=df.merge(cs[["protein","start","end","pen_nonhomolog"]],on=["protein","start","end"],how="left")
for m in df[df.pen_nonhomolog.isna()].itertuples():
    print("WARNING: no exact-match screen row for",m.protein,m.start,m.end,"- spec_other provisional (pen_nonhomolog=0)")
df["pen_nonhomolog"]=df.pen_nonhomolog.fillna(0)
pr=pd.read_csv("per_residue_calls.csv")
acc,tools,pl=[],[],[]
for r in df.itertuples():
    s=pr[(pr.protein==r.protein)&pr.pos.between(r.start,r.end)]
    acc.append(float((s.rsa.astype(float)>=0.25).mean())); tools.append(round(s.n_tools.mean(),2))
    pl.append(None if r.protein=="hlyE" else round(s.plddt.astype(float).mean(),1))
df["accessibility"]=acc; df["tools_per_res"]=tools; df["mean_plddt"]=pl
df["spec_nts"]=1-df.p_share_nts
df["spec_other"]=(1-df.p_share_conf)*(1-df.pen_nonhomolog)
crit=["spec_nts","spec_other","accessibility","conservation"]
if os.path.exists("ranking_inputs_extra.csv"):
    ex=pd.read_csv("ranking_inputs_extra.csv")[["protein","prior_evidence"]].drop_duplicates()
    df=df.merge(ex,on="protein",how="left")
    if df.prior_evidence.notna().all(): crit.append("prior_evidence")
if os.path.exists("prior_segments.csv") and "prior_evidence" in df.columns:
    ps=pd.read_csv("prior_segments.csv")
    def seg_prior(r):
        if r.protein not in set(ps.protein): return r.prior_evidence
        o=ps[(ps.protein==r.protein)&(ps.start<=r.end)&(ps.end>=r.start)]
        return float(o.prior_evidence.max()) if len(o) else 0.0
    df["prior_evidence"]=[seg_prior(r) for r in df.itertuples()]
w=np.array([W[c] for c in crit]); w/=w.sum(); X=df[crit].to_numpy(float)
df["score"]=(X*w).sum(1); df["rank"]=df.score.rank(ascending=False,method="min").astype(int)
cn=[c for c in crit if c!="prior_evidence"]; wn=np.array([W[c] for c in cn]); wn/=wn.sum()
df["score_no_prior"]=(df[cn].to_numpy(float)*wn).sum(1)
df["rank_no_prior"]=df.score_no_prior.rank(ascending=False,method="min").astype(int)
def sens(pm,n=1000):
    same=np.zeros(len(df)); top5=np.zeros(len(df))
    for _ in range(n):
        ww=np.array([W[c] for c in crit])*rng.uniform(1-pm,1+pm,len(crit)); ww/=ww.sum()
        rk=pd.Series((X*ww).sum(1)).rank(ascending=False,method="min").to_numpy()
        same+=(rk==df["rank"].to_numpy()); top5+=(rk<=5)
    return (100*same/n).round(0),(100*top5/n).round(0)
df["keeps_rank_pm10_pct"],df["in_top5_pm10_pct"]=sens(.10); df["keeps_rank_pm20_pct"],df["in_top5_pm20_pct"]=sens(.20)
df=df.sort_values(["rank","protein","start"]); df.round(3).to_csv("epitope_ranking_v3.csv",index=False)
print("criteria:",crit,"weights:",dict(zip(crit,w.round(3))))
print(df[["rank","protein","start","end","length","spec_nts","spec_other","accessibility","conservation","score","klass","prior_evidence","rank_no_prior","mean_plddt","keeps_rank_pm20_pct","in_top5_pm20_pct"]].round(2).to_string(index=False))
# consistency: carrier counts vs gene-level percentages
for p,f in ALN.items():
    a=read_aln(f); c=collections.Counter(GRP[n.split("_")[0]] for n in a)
    print(p,"unique carrier genomes:",dict(c),"| gene-level %:",gl.loc[p].to_dict())
