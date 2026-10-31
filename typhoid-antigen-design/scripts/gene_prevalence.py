"""Gene prevalence per genome group, counted on UNIQUE genomes.

usage (from the repository root): python3 -I scripts/gene_prevalence.py

Reads results/genomes/presence_long.tsv and writes results/genomes/gene_prevalence_unique.csv.
presence_long.tsv has duplicate rows for six genomes (the four S. Typhi genomes that were also
downloaded into the non-typhoidal folder, and two Enteritidis genomes). Counting rows instead of
genomes gives 133 NTS entries; the unique count is 131, so NTS hlyE/cdtB prevalence is 7.6%/8.4%
here versus 7.5%/8.3% in the original gene_level_summary.csv.
"""
import csv, collections
GROUP = {"typhi": "typhoidal", "paratyphiA": "typhoidal", "enteritidis": "NTS", "typhimurium": "NTS",
         "enterica": "NTS", "ecoli": "confounders", "shigella": "confounders",
         "klebsiella": "confounders", "pseudomonas": "confounders"}
genomes, carriers = collections.defaultdict(set), collections.defaultdict(set)
for r in csv.DictReader(open("results/genomes/presence_long.tsv"), delimiter="\t"):
    g = GROUP[r["group"]]
    genomes[g].add(r["accession"])
    if r["protein"]:
        carriers[(r["gene"], g)].add(r["accession"])
genes = sorted({k[0] for k in carriers}, key=["hlyE", "yncE", "cdtB", "fliC", "ompA", "ompC", "ompF"].index)
with open("results/genomes/gene_prevalence_unique.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["gene", "group", "n_genomes", "n_carriers", "pct"])
    for gene in genes:
        for g in ("typhoidal", "NTS", "confounders"):
            n, c = len(genomes[g]), len(carriers[(gene, g)])
            w.writerow([gene, g, n, c, round(100 * c / n, 1)])
print({g: len(v) for g, v in genomes.items()})
