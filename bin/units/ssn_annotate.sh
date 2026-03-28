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

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <annotation.csv> <taxonomy.tsv>"
    exit 1
fi

ANNOT_CSV="$1"
OUT_ANNOT="$2"

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
for COL in locus_tag class order; do
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

awk -F',' '
BEGIN {
    OFS="\t"
}
NR==1 {
    # Locate required columns by name
    for (i = 1; i <= NF; i++) {
        # Strip surrounding whitespace/quotes from header tokens
        h = $i; gsub(/^[[:space:]"]+|[[:space:]"]+$/, "", h)
        if (h == "locus_tag") l = i
        if (h == "class")     c = i
        if (h == "order")     o = i
    }
    if (!l || !c || !o) {
        print "[ERROR] Required columns (locus_tag, class, order) not found in CSV" \
            > "/dev/stderr"
        exit 1
    }
    print "id", "class", "order"
    next
}
{
    id    = $l
    class = $c
    order = $o

    # Strip surrounding whitespace/quotes from values
    gsub(/^[[:space:]"]+|[[:space:]"]+$/, "", id)
    gsub(/^[[:space:]"]+|[[:space:]"]+$/, "", class)
    gsub(/^[[:space:]"]+|[[:space:]"]+$/, "", order)

    # Skip rows with no meaningful id
    if (id == "" || id == "NA") next

    if (class == "") class = "Unknown"
    if (order == "") order = "Unknown"

    print id, class, order
}
' "$ANNOT_CSV" \
| sort -k1,1 -u \
> "$OUT_ANNOT"

ANNOT_COUNT=$(( $(wc -l < "$OUT_ANNOT") - 1 ))
log "Annotation entries written: $ANNOT_COUNT"

[ "$ANNOT_COUNT" -gt 0 ] || {
    log "WARNING: No annotation entries were written. Check that locus_tag values are not all empty/NA."
}

log "Annotation completed successfully -> $OUT_ANNOT"