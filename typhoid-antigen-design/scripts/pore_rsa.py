import copy
import pandas as pd
from Bio.PDB import PDBParser, ShrakeRupley
MAX={"ALA":106,"ARG":248,"ASN":157,"ASP":163,"CYS":135,"GLN":198,"GLU":194,"GLY":84,"HIS":184,"ILE":169,"LEU":164,"LYS":205,"MET":188,"PHE":197,"PRO":136,"SER":130,"THR":142,"TRP":227,"TYR":222,"VAL":142}
ONE={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
s=PDBParser(QUIET=True).get_structure("p","2wcd.pdb")
def rsa(keep,target):
    m=copy.deepcopy(s)[0]
    for c in list(m):
        if c.id not in keep: m.detach_child(c.id)
        else:
            for r in list(c):
                if r.id[0]!=" ": c.detach_child(r.id)
    ShrakeRupley().compute(m,level="R")
    return {r.id[1]:(r.resname,r.sasa/MAX[r.resname]) for r in m[target]}
PA=set("ABCDEFGHIJKL"); PM=set("MNOPQRSTUVWX")
R={"mono_A":rsa({"A"},"A"),"pore_A":rsa(PA,"A"),"mono_M":rsa({"M"},"M"),"pore_M":rsa(PM,"M")}
pos=sorted(R["mono_A"])
df=pd.DataFrame({"pos":pos})
df["aa"]=[ONE[R["mono_A"][p][0]] for p in pos]
for k,v in R.items(): df[k]=[round(v[p][1],3) for p in pos]
pr=pd.read_csv("per_residue_calls.csv"); pr=pr[pr.protein=="hlyE"]
df["rsa_1QOY"]=df.pos.map(pd.to_numeric(pr.set_index("pos").rsa,errors="coerce")).round(3)
df.to_csv("pore_rsa_residues.csv",index=False)
ce=pd.read_csv("consensus_epitopes.csv"); ce=ce[ce.protein=="hlyE"]
segs=sorted(set(zip(ce.start,ce.end))|{(54,61),(180,188)})
cols=["rsa_1QOY","mono_A","pore_A","mono_M","pore_M"]
rows=[]
for a,b in segs:
    d=df[df.pos.between(a,b)]
    rows.append(dict(start=a,end=b,n_res=len(d),**{c:round(float((d[c]>=0.25).mean()),2) for c in cols}))
out=pd.DataFrame(rows); out.to_csv("pore_rsa_segments.csv",index=False)
print("fraction of residues with RSA>=0.25"); print(out.to_string(index=False))
print(); print(df[df.pos.between(175,188)][["pos","aa"]+cols].to_string(index=False))
