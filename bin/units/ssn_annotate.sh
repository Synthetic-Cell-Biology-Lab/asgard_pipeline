#!/usr/bin/env bash
# =============================================================================
# SSN Step 5: Annotation — extract class & order from input CSV
# =============================================================================
# Usage:
#   ./ssn_annotate.sh <annotation.csv> <taxonomy.tsv>
#
# Arguments:
#   annotation.csv  Input CSV with at minimum these named columns:
#                     locus_tag  — sequence identifier matching FASTA headers
#                     class      — taxonomic or functional class
#                     order      — taxonomic or functional order
#   taxonomy.tsv    Output: TSV with columns (id, class, order)
#                   One row per unique locus_tag; blank values filled as "Unknown".
#                   Rows where locus_tag is empty or "NA" are skipped.
#
# Notes:
#   - Column order in the CSV does not matter; columns are found by header name.
#   - Duplicate locus_tag rows are deduplicated (first occurrence wins via sort -u).
#   - The output is sorted by id for stable diffs when parameters change.
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <annotation.csv> <taxonomy.tsv> <locus_tag_column_header>"
    exit 1
fi

ANNOT_CSV="$1"
OUT_ANNOT="$2"
ID_HEADER="$3"

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
die() { echo "[ERROR] $*" >&2; exit 1; }

check_tool() {
    command -v "$1" &>/dev/null || die "'$1' not found in PATH."
}

# -----------------------------------------------------------------------------
# Dependency checks
# -----------------------------------------------------------------------------

check_tool awk
check_tool sort

# -----------------------------------------------------------------------------
# Input validation
# -----------------------------------------------------------------------------

[ -f "$ANNOT_CSV" ] || die "Annotation CSV not found: $ANNOT_CSV"
[ -s "$ANNOT_CSV" ] || die "Annotation CSV is empty: $ANNOT_CSV"

# Quick sanity check: required column headers must be present in the first line
HEADER=$(head -1 "$ANNOT_CSV")
for COL in "$ID_HEADER" domain phylum class order; do
    echo "$HEADER" | grep -q "$COL" \
        || die "Required column '$COL' not found in CSV header: $ANNOT_CSV"
done

mkdir -p "$(dirname "$OUT_ANNOT")"

# -----------------------------------------------------------------------------
# Extract annotation columns dynamically
#
# Column positions are resolved from the header row so the CSV column order
# does not matter.  Rows with empty or "NA" locus_tag are skipped.
# Missing class/order values are filled with "Unknown".
# Duplicate ids are deduplicated via sort -k1,1 -u (first occurrence wins).
# -----------------------------------------------------------------------------

log "Extracting annotation (locus_tag, class, order) from $ANNOT_CSV..."
awk -F',' -v id_col="$ID_HEADER" '
BEGIN { OFS="\t" }

NR==1 {
    for (i = 1; i <= NF; i++) {
        h = $i
        gsub(/^[[:space:]"]+|[[:space:]"]+$/, "", h)

        if (h == id_col)   l = i
        if (h == "class")  c = i
        if (h == "order")  o = i
        if (h == "domain") d = i
        if (h == "phylum") p = i
    }

    if (!l || !c || !o || !d || !p) {
        print "[ERROR] Required columns missing" > "/dev/stderr"
        exit 1
    }

    header = "id\tdomain\tphylum\tclass\torder"
    print header
    next
}

{
    id     = $l
    domain = $d
    phylum = $p
    class  = $c
    order  = $o

    gsub(/^[[:space:]"]+|[[:space:]"]+$/, "", id)
    gsub(/^[[:space:]"]+|[[:space:]"]+$/, "", domain)
    gsub(/^[[:space:]"]+|[[:space:]"]+$/, "", phylum)
    gsub(/^[[:space:]"]+|[[:space:]"]+$/, "", class)
    gsub(/^[[:space:]"]+|[[:space:]"]+$/, "", order)

    if (id == "" || id == "NA") next

    if (domain == "") domain = "Unknown"
    if (phylum == "") phylum = "Unknown"
    if (class  == "") class  = "Unknown"
    if (order  == "") order  = "Unknown"

    print id, domain, phylum, class, order
}
' "$ANNOT_CSV" | {
    read -r header
    echo "$header"
    sort -k1,1 -u
} > "$OUT_ANNOT"

ANNOT_COUNT=$(( $(wc -l < "$OUT_ANNOT") - 1 ))
log "Annotation entries written: $ANNOT_COUNT"

[ "$ANNOT_COUNT" -gt 0 ] || {
    log "WARNING: No annotation entries were written. Check that locus_tag values are not all empty/NA."
}

log "Annotation completed successfully -> $OUT_ANNOT"