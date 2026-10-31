#!/usr/bin/env bash
# Run the ranking, classification and pore-accessibility scripts from the repository layout.
#
# The three scripts (rank_epitopes_v3.py, classify_segments.py, pore_rsa.py) read and write flat file
# names in the current directory, exactly as they were run during the analysis. They are kept unchanged;
# this wrapper stages the inputs into work/ (git-ignored) under those names, runs them there with
# python3 -I, and copies the outputs back into results/epitopes/.
#
# usage (from the repository root):  bash scripts/run_ranking_pipeline.sh [rank|classify|pore|all]   (ONE step per call; default: all)
# pore step needs structures/2wcd.pdb (PDB 2WCD, PDB format, downloaded from RCSB; not stored in the repo)
# and Biopython. Rank/classify need numpy and pandas.
set -euo pipefail
STEP="${1:-all}"
E=results/epitopes; G=results/genomes; D=data/input
rm -rf work && mkdir work

need() {   # need <source> <name expected by the scripts>
  if [ -f "$1" ]; then cp "$1" "work/$2"; else echo "MISSING: $1  (needed as $2)"; MISSING=1; fi
}
MISSING=0
if [ "$STEP" = rank ] || [ "$STEP" = classify ] || [ "$STEP" = all ]; then
  need $G/hlyE_aligned.aln aligned_hlye.aln
  need $G/yncE_aligned.aln yncE_aligned.aln
  need $G/cdtB_aligned.aln cdtB_aligned.aln
  need $D/full_length.faa full_length.faa
  need $G/gene_level_summary.csv gene_level_summary.csv
  need $E/consensus_epitopes.csv consensus_epitopes.csv
  need $E/per_residue_calls.csv per_residue_calls.csv
  need $E/epitope_crossreactivity_summary.csv epitope_crossreactivity_summary.csv
  need $D/prior_evidence_protein_level.csv ranking_inputs_extra.csv
  need $D/prior_segments.csv prior_segments.csv
  if [ "$STEP" = classify ]; then need $E/epitope_ranking_v3.csv epitope_ranking_v3.csv; fi   # classify reads the ranking table
fi
if [ "$STEP" = pore ] || [ "$STEP" = all ]; then
  need structures/2wcd.pdb 2wcd.pdb
  need $E/per_residue_calls.csv per_residue_calls.csv
  need $E/consensus_epitopes.csv consensus_epitopes.csv
fi
[ "$MISSING" = 0 ] || { echo "Add the missing files (see scripts/README.md) and rerun."; exit 1; }

cd work
if [ "$STEP" = rank ] || [ "$STEP" = all ]; then
  python3 -I ../scripts/rank_epitopes_v3.py > rank_run.txt && cp epitope_ranking_v3.csv ../$E/ && echo "ranking -> $E/epitope_ranking_v3.csv"
fi
if [ "$STEP" = classify ] || [ "$STEP" = all ]; then
  python3 -I ../scripts/classify_segments.py > classify_run.txt && cp epitope_classes.csv ../$E/ && echo "classes -> $E/epitope_classes.csv"
fi
if [ "$STEP" = pore ] || [ "$STEP" = all ]; then
  python3 -I ../scripts/pore_rsa.py > pore_run.txt && cp pore_rsa_segments.csv pore_rsa_residues.csv ../$E/ && echo "pore -> $E/pore_rsa_segments.csv, pore_rsa_residues.csv"
fi
