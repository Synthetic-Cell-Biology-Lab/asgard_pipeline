#!/usr/bin/env bash
# =============================================================================
# SSN Step 3: Filter edges + extract node metadata
# =============================================================================
# Usage:
#   ./ssn_filter.sh <similarities.tsv> <nr.fasta> <edges.tsv> <nodes.tsv> \
#                   <bitscore> <coverage> <evalue>
#
# Arguments:
#   similarities.tsv  Raw all-vs-all results (output of ssn_search.sh)
#   nr.fasta          Non-redundant FASTA (output of ssn_cdhit.sh)
#   edges.tsv         Output: filtered edge list (node1, node2, bitscore)
#   nodes.tsv         Output: node metadata table (id, full_header, length)
#   bitscore          Minimum bitscore to retain an edge (e.g. 50)
#   coverage          Minimum query/subject coverage fraction (e.g. 0.6)
#   evalue            Maximum e-value to retain an edge (e.g. 1e-5)
#
# Search engine is auto-detected from the similarities file:
#   BLAST output   has qlen/slen in cols 7-8  (raw alignment lengths)
#   MMseqs2 output has qcov/tcov in cols 7-8  (pre-computed fractions)
#
# The heuristic: if the median value of col7 > 1 it must be a raw length
# (fractions are always <= 1), so we treat it as BLAST output.
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <similarities.tsv> <edges.tsv> <bitscore> <coverage> <evalue>"
    exit 1
fi

SIMILARITIES="$1"
EDGES="$2"
BITSCORE="$3"
COVERAGE="$4"
EVALUE="$5"
NODES="$6"

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

[ -f "$SIMILARITIES" ] || die "Similarities file not found: $SIMILARITIES"
[ -s "$SIMILARITIES" ] || die "Similarities file is empty: $SIMILARITIES"


for F in "$EDGES"; do
    mkdir -p "$(dirname "$F")"
done

# -----------------------------------------------------------------------------
# Auto-detect search method from col7 values
# Fractions (MMseqs2) are always <= 1; raw lengths (BLAST) are > 1
# Sample first 200 data lines for the check
# -----------------------------------------------------------------------------

METHOD=$(awk '
    NR <= 200 {
        val = $7 + 0
        if (val > 1) { blast++; next }
        mmseqs++
    }
    END {
        if (blast > mmseqs) print "blast"
        else                print "mmseqs"
    }
' "$SIMILARITIES")

log "Detected search method: ${METHOD}"
log "Filtering edges (bitscore>=${BITSCORE}, coverage>=${COVERAGE}, evalue<=${EVALUE})..."

# -----------------------------------------------------------------------------
# Filter edges
#
# BLAST columns:   qseqid sseqid pident alnlen evalue bitscore qlen   slen
#                  1      2      3      4      5      6        7      8
#   coverage = alnlen / qlen  and  alnlen / slen
#
# MMseqs2 columns: query  target pident alnlen evalue bits     qcov   tcov
#                  1      2      3      4      5      6        7      8
#   coverage already provided as fractions in cols 7 & 8
#
# In both cases self-hits ($1 == $2) are discarded.
# Edges are deduplicated by always printing the lexicographically smaller
# node first and then running sort -u.
# -----------------------------------------------------------------------------

if [ "$METHOD" = "blast" ]; then
    awk \
        -v bits="$BITSCORE"   \
        -v cov="$COVERAGE"    \
        -v max_eval="$EVALUE" \
        'BEGIN { OFS="\t" }
        {
            if ($1 == $2)              next
            if ($6 < bits)             next
            if ($5 + 0 > max_eval + 0) next
            qcov = ($4 + 0) / ($7 + 0)
            scov = ($4 + 0) / ($8 + 0)
            if (qcov < cov || scov < cov) next
            if ($1 < $2) print $1, $2, $6
            else         print $2, $1, $6
        }' "$SIMILARITIES" \
    | sort -k1,1 -k2,2 -u \
    > "$EDGES"
else
    awk \
        -v bits="$BITSCORE"   \
        -v cov="$COVERAGE"    \
        -v max_eval="$EVALUE" \
        'BEGIN { OFS="\t" }
        {
            if ($1 == $2)              next
            if ($6 < bits)             next
            if ($5 + 0 > max_eval + 0) next
            if ($7 < cov || $8 < cov)  next
            if ($1 < $2) print $1, $2, $6
            else         print $2, $1, $6
        }' "$SIMILARITIES" \
    | sort -k1,1 -k2,2 -u \
    > "$EDGES"
fi

EDGE_COUNT=$(wc -l < "$EDGES")
log "Edges after filtering: $EDGE_COUNT"

if [ "$EDGE_COUNT" -eq 0 ]; then
    log "WARNING: No edges passed filters."
    log "  Consider relaxing: bitscore (current: ${BITSCORE}), coverage (current: ${COVERAGE}), evalue (current: ${EVALUE})"
fi

NODE_COUNT=$(( $(wc -l < "$NODES") - 1 ))
log "Nodes in metadata table: $NODE_COUNT"


# -----------------------------------------------------------------------------
# Network statistics
# -----------------------------------------------------------------------------

log "Network statistics:"

if [ "$EDGE_COUNT" -gt 0 ]; then
    CONNECTED_NODES=$(awk 'NR>1{print $1"\n"$2}' "$EDGES" | sort -u | wc -l)
    ISOLATED_NODES=$(( NODE_COUNT - CONNECTED_NODES ))
    AVG_DEGREE=$(awk -v e="$EDGE_COUNT" -v n="$CONNECTED_NODES" \
                     'BEGIN { printf "%.2f", (2*e)/n }')
    log "  Connected nodes : $CONNECTED_NODES / $NODE_COUNT"
    log "  Isolated nodes  : $ISOLATED_NODES"
    log "  Edges           : $EDGE_COUNT"
    log "  Mean degree     : $AVG_DEGREE"
else
    log "  No edges -- network is empty."
fi

log "Filter completed successfully."
log "  Edges -> $EDGES"
