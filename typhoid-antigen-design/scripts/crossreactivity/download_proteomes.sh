#!/usr/bin/env bash
# Download UniProt proteomes listed in proteomes.tsv -> data/proteomes/<name>.faa
# usage: bash download_proteomes.sh [proteomes.tsv] [outdir]
set -euo pipefail
TSV=${1:-proteomes.tsv}; OUT=${2:-data/proteomes}; mkdir -p "$OUT"
while IFS=$'\t' read -r name taxid upid group; do
  [[ -z "${name:-}" || "$name" == \#* ]] && continue
  if [[ -z "${upid:-}" ]]; then
    upid=$(curl -fsS "https://rest.uniprot.org/proteomes/search?query=(organism_id:${taxid})%20AND%20(proteome_type:1)&fields=upid&format=tsv&size=1" | sed -n '2p' || true)
    [[ -z "$upid" ]] && { echo "!! $name: no reference proteome found for taxid $taxid - put a UP... ID in the TSV" >&2; continue; }
  fi
  echo ">> $name  taxid=$taxid  proteome=$upid  group=$group"
  curl -fsS "https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=(proteome:${upid})" -o "$OUT/$name.faa"
  n=$(grep -c '^>' "$OUT/$name.faa"); echo "   $n proteins; first header: $(head -1 "$OUT/$name.faa" | cut -c1-110)"
done < "$TSV"
echo "Check each first header above matches the organism you intended."
