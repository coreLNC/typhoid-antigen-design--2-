#!/usr/bin/env bash
# Whole-protein DIAMOND screen of each antigen against each proteome.
# usage: bash run_diamond.sh full_length.faa [data/proteomes] [results/diamond]   (run inside: conda activate antigen)
set -euo pipefail
Q=${1:-full_length.faa}; P=${2:-data/proteomes}; O=${3:-results/diamond}; mkdir -p "$O" "$P/db"
for f in "$P"/*.faa; do
  n=$(basename "$f" .faa)
  [[ -f "$P/db/$n.dmnd" ]] || diamond makedb --in "$f" -d "$P/db/$n" --quiet
  diamond blastp -q "$Q" -d "$P/db/$n" --evalue 1e-3 --more-sensitive --max-target-seqs 5 \
    --outfmt 6 qseqid sseqid pident length qcovhsp evalue bitscore stitle -o "$O/$n.tsv" --quiet
  echo "$n: $(wc -l < "$O/$n.tsv") hits"
done
