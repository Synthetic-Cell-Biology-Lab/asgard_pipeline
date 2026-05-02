#!/usr/bin/env bash

set -euo pipefail

CSV="$1"          # taxonomy CSV
FASTA="$2"        # protein fasta
LEVEL="$3"        # taxonomic level (e.g. genus, family, phylum)
OUTDIR="${4:-msa_by_taxon}"

mkdir -p "$OUTDIR/tmp"

echo "[INFO] Parsing CSV and grouping by $LEVEL..."

# Get column index of the requested level
HEADER=$(head -n1 "$CSV")
IFS=',' read -ra COLS <<< "$HEADER"

LEVEL_IDX=-1
LOCUS_IDX=-1

for i in "${!COLS[@]}"; do
    col=$(echo "${COLS[$i]}" | tr -d '"')
    if [[ "$col" == "$LEVEL" ]]; then
        LEVEL_IDX=$((i+1))
    fi
    if [[ "$col" == "locus_tag" ]]; then
        LOCUS_IDX=$((i+1))
    fi
done

if [[ "$LEVEL_IDX" -lt 0 || "$LOCUS_IDX" -lt 0 ]]; then
    echo "[ERROR] Could not find required columns"
    exit 1
fi

# Build mapping: taxon → locus_tags
awk -F',' -v lvl="$LEVEL_IDX" -v loc="$LOCUS_IDX" '
NR>1 {
    gsub(/"/,"",$lvl);
    gsub(/"/,"",$loc);
    if($lvl != "" && $loc != "") {
        print $lvl "\t" $loc
    }
}' "$CSV" | sort | uniq > "$OUTDIR/tmp/taxon_map.tsv"

echo "[INFO] Creating per-taxon FASTA files..."

# Loop over taxa
cut -f1 "$OUTDIR/tmp/taxon_map.tsv" | sort | uniq | while read -r TAXON; do
    SAFE_TAXON=$(echo "$TAXON" | tr ' /' '__')
    LIST="$OUTDIR/tmp/${SAFE_TAXON}.ids"
    OUTFA="$OUTDIR/${SAFE_TAXON}.fasta"

    # Extract locus_tags for this taxon
    awk -v t="$TAXON" '$1==t {print $2}' "$OUTDIR/tmp/taxon_map.tsv" > "$LIST"

    # Subset FASTA
    awk '
    BEGIN {
        while((getline line < "'$LIST'") > 0) {
            ids[line]=1
        }
    }
    /^>/ {
        header=$0
        id=$1
        sub(/^>/,"",id)
        keep=ids[id]
    }
    keep { print }
    ' "$FASTA" > "$OUTFA"

    # Skip tiny groups
    SEQ_COUNT=$(grep -c "^>" "$OUTFA" || true)
    if [[ "$SEQ_COUNT" -lt 2 ]]; then
        rm -f "$OUTFA"
        continue
    fi

    echo "[INFO] Running MSA for $TAXON ($SEQ_COUNT seqs)..."

    # Run MAFFT (you can swap with Clustal Omega)
    mafft --auto "$OUTFA" > "$OUTDIR/${SAFE_TAXON}.aln.fasta"

done

echo "[DONE] MSAs written to $OUTDIR"