#!/usr/bin/env python3
"""Epitope-level cross-reactivity screen (exact k-mer / longest shared substring) + DIAMOND summary.
usage: python3 -I screen_epitopes.py consensus_epitopes.csv proteomes.tsv data/proteomes results [diamond_dir]
Outputs: results/epitope_crossreactivity_long.csv, results/epitope_crossreactivity_summary.csv,
         results/antigen_diamond_summary.csv (if diamond_dir given)"""
import sys, os, bisect, glob
import pandas as pd
KS=(5,6,7,8); RISK_LEN=8     # exact shared stretch >= RISK_LEN residues (or the whole segment, if shorter) flags HIGH risk
def read_fasta(p):
    ids,descs,seqs,cur=[],[],[],[]
    for line in open(p):
        line=line.rstrip()
        if line.startswith(">"):
            if ids: seqs.append("".join(cur)); cur=[]
            h=line[1:]; ids.append(h.split()[0]); descs.append(h)
        else: cur.append(line)
    if ids: seqs.append("".join(cur))
    return ids,descs,seqs
class Proteome:
    def __init__(s,path):
        s.ids,s.descs,seqs=read_fasta(path)
        s.text="*".join(seqs).upper(); s.starts=[]; o=0
        for q in seqs: s.starts.append(o); o+=len(q)+1
    def proteins_with(s,kmer):
        found=set(); i=s.text.find(kmer)
        while i!=-1:
            found.add(bisect.bisect_right(s.starts,i)-1); i=s.text.find(kmer,i+1)
        return found
def screen(pep,prot):
    pep=pep.upper(); out={}
    for k in KS:
        wins=[pep[i:i+k] for i in range(len(pep)-k+1)]
        out[f"frac_{k}mer_hit"]=round(sum(bool(prot.proteins_with(w)) for w in wins)/len(wins),3) if wins else None
    best=(0,"",None)
    for L in range(5,len(pep)+1):
        hit=None
        for i in range(len(pep)-L+1):
            ps=prot.proteins_with(pep[i:i+L])
            if ps: hit=(L,pep[i:i+L],min(ps)); break
        if hit is None: break
        best=hit
    out["longest_exact_match"]=best[0] if best[0] else 4   # <5 means no 5-mer shared
    out["longest_match_peptide"]=best[1]
    out["best_protein"]=prot.descs[best[2]][:90] if best[2] is not None else ""
    return out
def main():
    epi,tsv,pdir,outdir=sys.argv[1:5]; ddir=sys.argv[5] if len(sys.argv)>5 else None
    os.makedirs(outdir,exist_ok=True)
    epi=pd.read_csv(epi)
    cfg=[l.rstrip("\n").split("\t") for l in open(tsv) if l.strip() and not l.startswith("#")]
    rows=[]
    for name,taxid,upid,group in [(c[0],c[1],c[2],c[3]) for c in cfg]:
        f=os.path.join(pdir,name+".faa")
        if not os.path.exists(f): print("skip (not downloaded):",name); continue
        prot=Proteome(f); print(f"{name}: {len(prot.ids)} proteins, {len(prot.text):,} residues")
        for _,e in epi.iterrows():
            r=screen(e.peptide,prot); r.update(protein=e.protein,start=e.start,end=e.end,length=e.length,peptide=e.peptide,proteome=name,group=group)
            rows.append(r)
    long=pd.DataFrame(rows); first=["protein","start","end","length","peptide","proteome","group"]
    long=long[first+[c for c in long.columns if c not in first]]
    long.to_csv(os.path.join(outdir,"epitope_crossreactivity_long.csv"),index=False)
    conf=long[long.group!="typhoidal_related"]
    summ=[]
    for (p,a,b),g in conf.groupby(["protein","start","end"],sort=False):
        w=g.loc[g.longest_exact_match.idxmax()]
        row=dict(protein=p,start=a,end=b,length=int(g.length.iloc[0]),peptide=g.peptide.iloc[0],
                 worst_proteome=w.proteome,worst_group=w.group,longest_exact_match=int(w.longest_exact_match),
                 shared_peptide=w.longest_match_peptide,worst_hit=w.best_protein,
                 max_frac_7mer_hit=g.frac_7mer_hit.max(),
                 risk=("HIGH" if (w.longest_exact_match>=RISK_LEN or (len(g.peptide.iloc[0])>=7 and w.longest_exact_match>=len(g.peptide.iloc[0]))) else "low"))
        for grp,gg in g.groupby("group"): row[f"longest_{grp}"]=int(gg.longest_exact_match.max())
        summ.append(row)
    summ=pd.DataFrame(summ); summ.to_csv(os.path.join(outdir,"epitope_crossreactivity_summary.csv"),index=False)
    print("\n",summ.drop(columns=["peptide","worst_hit"]).to_string(index=False))
    if ddir:
        d=[]
        for f in glob.glob(os.path.join(ddir,"*.tsv")):
            if os.path.getsize(f)==0: continue
            t=pd.read_csv(f,sep="\t",header=None,names=["antigen","hit","pident","aln_len","qcov","evalue","bitscore","stitle"])
            t["proteome"]=os.path.basename(f)[:-4]; d.append(t.sort_values("bitscore",ascending=False).groupby("antigen").head(1))
        if d:
            d=pd.concat(d); d.to_csv(os.path.join(outdir,"antigen_diamond_summary.csv"),index=False)
            print("\nDIAMOND best hit per antigen/proteome:\n",d[["antigen","proteome","pident","qcov","evalue","stitle"]].to_string(index=False))
if __name__=="__main__": main()
