# Scripts

| Script | Purpose | Status |
|---|---|---|
| `gene_prevalence.py` | gene prevalence per group on unique genomes | here |
| `make_figures.py` | Figures 1, 4, 6, 7, 8 | here |
| `render_structures.py` | Figure 5 (3Dmol.js renders) | here |
| `extract_west_african_genomes.py` | West African genome table from NCBI Datasets JSON | here |
| `consensus.py` | three-tool consensus epitopes (BepiPred-3.0, DiscoTope-3.0, ElliPro) | here |
| `crossreactivity/` | DIAMOND and exact k-mer screens against NTS, confounders, human, *Plasmodium* | here |
| `archive/rank_epitopes_v2.py` | earlier ranking script (superseded by v3) | here, historical |
| `rank_epitopes_v3.py` | current ranking script (primary no-prior and prior-weighted ranks, weight perturbation); unchanged from the analysis run (its docstring still says v2) | here |
| `classify_segments.py` | pan-typhoidal vs Typhi-selective classification | here |
| `pore_rsa.py` | monomer vs pore accessibility (1QOY, 2WCD) | here; output verified identical to `pore_rsa_segments.csv` |
| `run_ranking_pipeline.sh` | stages inputs under the flat names the three scripts expect, runs them in `work/`, copies results back | here |
| genome, presence, alignment and sliding-window scripts (`scripts/pipeline/`) | produce `presence_long.tsv`, the `*_aligned.aln` and `*_windows_by_group.csv` files | copy from `~/typhoid_antigen/` |

All scripts are run from the repository root and use relative paths. `python3 -I` is used so that
no file in the working directory can shadow a standard-library module.

## Inputs still to add for the ranking and classification steps

`run_ranking_pipeline.sh` lists any missing input. Copy these from the analysis machine:

- `results/genomes/yncE_aligned.aln`
- `results/epitopes/epitope_crossreactivity_summary.csv` (30 rows, includes HlyE 175-188)
- `data/input/prior_segments.csv` (segment-level prior evidence, from Chin et al. 2017)
- `structures/2wcd.pdb` for the pore step (PDB 2WCD in PDB format from RCSB; git-ignored because of its size)

`results/epitopes/consensus_epitopes.csv` carries the appended, flagged literature row `hlyE,175,188` with blank tool columns,
as used for the v3 ranking.
