"""Three-tool B-cell epitope consensus (BepiPred-3.0, DiscoTope-3.0, ElliPro), full-length numbering.
usage: python3 -I consensus.py     (reads bepipred_raw_output.csv, ellipro_raw/*.txt, discotope_raw/*.csv)"""
import re, pandas as pd
full={r.split("\n")[0].replace("_full",""):"".join(r.split("\n")[1:]) for r in open("full_length.faa").read().split(">")[1:]}
raw=pd.read_csv("bepipred_raw_output.csv",skipinitialspace=True); raw.columns=["acc","aa","score","lin"]
starts={"hlyE":1,"yncE":31,"cdtB":23}; THR=0.1512; MIN_TOOLS=2
def ellipro(p):
    t=open(f"ellipro_raw/{p}.txt").read()
    lin_part,disc_part=t.split("Predicted Discontinuous Epitope(s)")
    lin=set()
    for m in re.finditer(r"^\s*\d+\s+A\s+(\d+)\s+(\d+)\s+[A-Z]+\s+\d+",lin_part,re.M):
        lin|=set(range(int(m[1]),int(m[2])+1))
    return lin,set(int(x) for x in re.findall(r"A:[A-Z](\d+)",disc_part))
rows=[]; mism={}
for p in full:
    acc=[a for a in raw.acc.unique() if a.startswith(p)][0]
    g=raw[raw.acc==acc].reset_index(drop=True); g["pos"]=g.index+starts[p]
    lin,disc=ellipro(p)
    dt=pd.read_csv(f"discotope_raw/{p}.csv").set_index("res_id")
    for _,r in g.iterrows():
        pos=int(r.pos); d=dt.loc[pos] if pos in dt.index else None
        if d is not None and d.residue!=r.aa:
            assert p=='hlyE', (p,pos,d.residue,r.aa)   # YncE/CdtB must match exactly; HlyE structure is E. coli (~91% identical)
            mism[p]=mism.get(p,0)+1
        rows.append(dict(protein=p,pos=pos,aa=r.aa,bepipred=bool(r.score>=THR),bepi_score=round(r.score,3),
            ellipro=pos in lin or pos in disc, ellipro_linear=pos in lin, ellipro_disc=pos in disc,
            discotope=bool(d.epitope) if d is not None else False, disco_calibrated=None if d is None else d.calibrated_score,
            rsa=None if d is None else d.rsa, plddt=None if d is None or p=="hlyE" else d.pLDDTs))
print("residue mismatches vs structure:",mism or "none"); df=pd.DataFrame(rows); df["n_tools"]=df[["bepipred","ellipro","discotope"]].sum(axis=1)
df.to_csv("per_residue_calls.csv",index=False)
print(df.groupby("protein").agg(residues=("pos","size"),bepipred=("bepipred","sum"),ellipro=("ellipro","sum"),
      discotope=("discotope","sum"),ge2=("n_tools",lambda s:(s>=2).sum()),all3=("n_tools",lambda s:(s==3).sum())).to_string())
def segs(sub,col,gap=1,minlen=4):
    v=sub[col].tolist(); pos=sub.pos.tolist(); out=[]; i=0
    while i<len(v):
        if v[i]:
            last=i;k=i+1;g=0
            while k<len(v):
                if v[k]: last=k;g=0
                else:
                    g+=1
                    if g>gap: break
                k+=1
            if last-i+1>=minlen: out.append((pos[i],pos[last]))
            i=last+1
        else: i+=1
    return out
df["call"]=df.n_tools>=MIN_TOOLS
res=[]
for p in full:
    sub=df[df.protein==p].reset_index(drop=True)
    for a,b in segs(sub,"call"):
        s=sub[(sub.pos>=a)&(sub.pos<=b)]
        res.append(dict(protein=p,start=a,end=b,length=b-a+1,peptide="".join(s.aa),
            mean_tools_per_res=round(s.n_tools.mean(),2),pct_all3=round(100*(s.n_tools==3).mean()),
            pct_bepipred=round(100*s.bepipred.mean()),pct_discotope=round(100*s.discotope.mean()),pct_ellipro=round(100*s.ellipro.mean()),
            mean_bepi=round(s.bepi_score.mean(),3),mean_rsa=round(s.rsa.astype(float).mean(),2),
            mean_plddt=None if p=="hlyE" else round(s.plddt.astype(float).mean(),1)))
out=pd.DataFrame(res); out.to_csv("consensus_epitopes.csv",index=False)
print(); print(out.drop(columns="peptide").to_string(index=False))
