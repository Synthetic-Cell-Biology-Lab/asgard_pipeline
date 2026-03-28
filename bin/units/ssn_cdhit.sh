#!/usr/bin/env bash
# =============================================================================
# SSN Step 1: CD-HIT clustering — remove redundant sequences
# =============================================================================
# Usage:
#   ./ssn_cdhit.sh <input.fasta> <nr.fasta> <identity> <threads>
#
# Arguments:
#   input.fasta   Raw input FASTA (protein sequences)
#   nr.fasta      Output: non-redundant FASTA
#   identity      Sequence identity threshold (e.g. 0.90)
#   threads       Number of CPU threads
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input.fasta> <nr.fasta> <identity> <threads>"
    exit 1
fi

INPUT="$1"
NR_FASTA="$2"
IDENTITY="$3"
THREADS="$4"

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

log "Checking dependencies..."
check_tool cd-hit
check_tool grep

# -----------------------------------------------------------------------------
# Input validation
# -----------------------------------------------------------------------------

[ -f "$INPUT" ]       || die "Input FASTA not found: $INPUT"
grep -q "^>" "$INPUT" || die "Input does not look like a FASTA file: $INPUT"

RAW_COUNT=$(grep -c "^>" "$INPUT")
log "Input sequences: $RAW_COUNT"
[ "$RAW_COUNT" -ge 2 ] || die "Need at least 2 sequences to build a network."

mkdir -p "$(dirname "$NR_FASTA")"

# -----------------------------------------------------------------------------
# CD-HIT
# -----------------------------------------------------------------------------

log "Running CD-HIT (identity=${IDENTITY}, threads=${THREADS})..."

cd-hit \
    -i "$INPUT"    \
    -o "$NR_FASTA" \
    -c "$IDENTITY" \
    -n 5           \
    -T "$THREADS"  \
    -M 0           \
    -d 0           \
    -g 1           \
    > "${NR_FASTA}.cdhit.log" 2>&1

# CD-HIT does not reliably exit non-zero on failure — check output explicitly
[ -f "$NR_FASTA" ]       || die "CD-HIT produced no output FASTA. See ${NR_FASTA}.cdhit.log"
[ -s "$NR_FASTA" ]       || die "CD-HIT output FASTA is empty. See ${NR_FASTA}.cdhit.log"
grep -q "^>" "$NR_FASTA" || die "CD-HIT output is not a valid FASTA. See ${NR_FASTA}.cdhit.log"

SEQ_COUNT=$(grep -c "^>" "$NR_FASTA")
log "Sequences after CD-HIT: $SEQ_COUNT (removed $((RAW_COUNT - SEQ_COUNT)))"
[ "$SEQ_COUNT" -ge 2 ] || die "Fewer than 2 sequences remain after CD-HIT. Lower the identity threshold."

log "CD-HIT completed successfully -> $NR_FASTA"