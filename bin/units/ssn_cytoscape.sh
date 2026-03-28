#!/usr/bin/env bash
# =============================================================================
# SSN Step 4: Cytoscape export — write SIF and edge-attribute files
# =============================================================================
# Usage:
#   ./ssn_cytoscape.sh <edges.tsv> <output.sif> <output.ea>
#
# Arguments:
#   edges.tsv     Filtered edge list (output of ssn_filter.sh)
#                 Columns: node1 <TAB> node2 <TAB> bitscore
#   output.sif    Output: Cytoscape Simple Interaction Format
#                 Columns: node1 <TAB> similarity <TAB> node2
#   output.ea     Output: Cytoscape edge-attribute file
#                 Contains the bitscore for each edge in "node1 (similarity) node2" format
#
# To load in Cytoscape:
#   1. File > Import > Network from File  -> <output.sif>
#   2. File > Import > Table from File    -> <output.ea>  (edge bitscore attribute)
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <edges.tsv> <output.sif> <output.ea>"
    exit 1
fi

EDGES="$1"
SIF="$2"
EA="$3"

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

# -----------------------------------------------------------------------------
# Input validation
# -----------------------------------------------------------------------------

[ -f "$EDGES" ] || die "Edges file not found: $EDGES"

for F in "$SIF" "$EA"; do
    mkdir -p "$(dirname "$F")"
done

EDGE_COUNT=$(wc -l < "$EDGES")
log "Edges to export: $EDGE_COUNT"

if [ "$EDGE_COUNT" -eq 0 ]; then
    log "WARNING: Edge file is empty. Cytoscape files will be written but will contain no edges."
fi

# -----------------------------------------------------------------------------
# SIF file
# Format: node1 <TAB> interaction_type <TAB> node2
# -----------------------------------------------------------------------------

log "Writing SIF -> $SIF"
awk 'BEGIN{OFS="\t"} {print $1, "similarity", $2}' "$EDGES" > "$SIF"

# -----------------------------------------------------------------------------
# Edge attribute file
# Header line: attribute name
# Data lines:  node1 (interaction_type) node2 = value
# -----------------------------------------------------------------------------

log "Writing edge attributes -> $EA"
printf "bitscore\n" > "$EA"
awk 'BEGIN{OFS="\t"} {print $1" (similarity) "$2, $3}' "$EDGES" >> "$EA"

log "Cytoscape export completed successfully."
log "  SIF -> $SIF"
log "  EA  -> $EA"