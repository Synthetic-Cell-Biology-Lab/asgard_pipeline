#!/usr/bin/env bash
# =============================================================================
# SSN Step 2: All-vs-all similarity search (BLAST or MMseqs2)
# =============================================================================
# Usage:
#   ./ssn_search.sh <nr.fasta> <similarities.tsv> <threads> <evalue>
#
# Arguments:
#   nr.fasta          Non-redundant FASTA (output of ssn_cdhit.sh)
#   similarities.tsv  Output: raw all-vs-all search results
#   threads           Number of CPU threads
#   evalue            Max e-value threshold (e.g. 1e-5)
#
# Search engine selection:
#   < 5000 sequences  -> BLAST
#   >= 5000 sequences -> MMseqs2
#
# Output columns (both methods):
#   1:query/qseqid  2:target/sseqid  3:pident  4:alnlen
#   5:evalue        6:bitscore       7:qcov     8:tcov
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <nr.fasta> <similarities.tsv> <threads> <evalue>"
    exit 1
fi

NR_FASTA="$1"
SIMILARITIES="$2"
THREADS="$3"
EVALUE="$4"

BLAST_THRESHOLD=5000

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
die() { echo "[ERROR] $*" >&2; exit 1; }

check_tool() {
    command -v "$1" &>/dev/null || die "'$1' not found in PATH."
}

# -----------------------------------------------------------------------------
# Input validation
# -----------------------------------------------------------------------------

[ -f "$NR_FASTA" ]       || die "NR FASTA not found: $NR_FASTA"
grep -q "^>" "$NR_FASTA" || die "NR FASTA does not look like a FASTA file: $NR_FASTA"

SEQ_COUNT=$(grep -c "^>" "$NR_FASTA")
log "Sequences to search: $SEQ_COUNT"
[ "$SEQ_COUNT" -ge 2 ] || die "Need at least 2 sequences to run all-vs-all search."

mkdir -p "$(dirname "$SIMILARITIES")"

# -----------------------------------------------------------------------------
# Select search engine
# -----------------------------------------------------------------------------

if [ "$SEQ_COUNT" -lt "$BLAST_THRESHOLD" ]; then
    METHOD="blast"
else
    METHOD="mmseqs"
fi

log "Search engine: ${METHOD} (threshold=${BLAST_THRESHOLD})"

# -----------------------------------------------------------------------------
# BLAST branch
# -----------------------------------------------------------------------------

if [ "$METHOD" = "blast" ]; then

    check_tool makeblastdb
    check_tool blastp

    log "Building BLAST database..."
    makeblastdb \
        -in     "$NR_FASTA"           \
        -dbtype prot                  \
        -out    "${NR_FASTA}.blastdb" \
        > /dev/null 2>&1

    log "Running blastp (evalue=${EVALUE}, threads=${THREADS})..."
    # Output columns:
    #   1:qseqid  2:sseqid  3:pident  4:length  5:evalue  6:bitscore  7:qlen  8:slen
    # Coverage is computed from alnlen/qlen and alnlen/slen in ssn_filter.sh
    blastp \
        -query       "$NR_FASTA"           \
        -db          "${NR_FASTA}.blastdb" \
        -out         "$SIMILARITIES"       \
        -outfmt      "6 qseqid sseqid pident length evalue bitscore qlen slen" \
        -num_threads "$THREADS"            \
        -evalue      "$EVALUE"             \
        -seg         yes

# -----------------------------------------------------------------------------
# MMseqs2 branch
# -----------------------------------------------------------------------------

else

    check_tool mmseqs

    TMPDIR=$(mktemp -d)
    trap 'rm -rf "$TMPDIR"' EXIT

    MMSEQS_WD="$TMPDIR/mmseqs_work"
    mkdir -p "$MMSEQS_WD"

    log "Creating MMseqs2 database..."
    mmseqs createdb "$NR_FASTA" "$TMPDIR/seqDB" > /dev/null 2>&1

    log "Running MMseqs2 search (sensitivity=7.5, evalue=${EVALUE}, threads=${THREADS})..."
    mmseqs search \
        "$TMPDIR/seqDB"    \
        "$TMPDIR/seqDB"    \
        "$TMPDIR/resultDB" \
        "$MMSEQS_WD"       \
        --threads        "$THREADS" \
        -s               7.5        \
        --num-iterations 2          \
        -e               "$EVALUE"  \
        > /dev/null 2>&1

    log "Converting MMseqs2 results..."
    # Output columns:
    #   1:query  2:target  3:pident  4:alnlen  5:evalue  6:bits  7:qcov  8:tcov
    # qcov/tcov are already fractions [0-1] — no conversion needed in ssn_filter.sh
    mmseqs convertalis \
        "$TMPDIR/seqDB"    \
        "$TMPDIR/seqDB"    \
        "$TMPDIR/resultDB" \
        "$SIMILARITIES"    \
        --format-output "query,target,pident,alnlen,evalue,bits,qcov,tcov" \
        > /dev/null 2>&1

fi

# -----------------------------------------------------------------------------
# Validate output
# -----------------------------------------------------------------------------

RAW_HITS=$(wc -l < "$SIMILARITIES")
log "Raw similarity hits: $RAW_HITS"
[ "$RAW_HITS" -gt 0 ] || die "No similarity hits found. Check your input or search parameters."

log "Search completed successfully -> $SIMILARITIES"