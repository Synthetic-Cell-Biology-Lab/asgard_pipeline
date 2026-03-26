#!/usr/bin/env bash
# =============================================================================
# Sequence Similarity Network (SSN) Pipeline
# =============================================================================
# Usage:
#   ./ssn_pipeline.sh <input.fasta> <nr.fasta> <similarities.tsv> \
#                     <edges.tsv> <nodes.tsv> <threads> <identity> \
#                     <bitscore> <coverage> [evalue]
#
# Arguments:
#   input.fasta      Raw input FASTA (protein sequences)
#   nr.fasta         Output: non-redundant FASTA after CD-HIT
#   similarities.tsv Output: raw all-vs-all search results
#   edges.tsv        Output: filtered edges (node1, node2, bitscore)
#   nodes.tsv        Output: node metadata table
#   threads          Number of CPU threads
#   identity         CD-HIT sequence identity threshold (e.g. 0.9)
#   bitscore         Minimum bitscore to retain an edge (e.g. 50)
#   coverage         Minimum query/subject coverage fraction (e.g. 0.8)
#   evalue           [optional] Max e-value threshold (default: 1e-5)
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------

if [ "$#" -lt 9 ]; then
    echo "Usage: $0 <input> <nr_fasta> <similarities> <edges> <nodes> <threads> <identity> <bitscore> <coverage> [evalue]"
    exit 1
fi

INPUT="$1"
NR_FASTA="$2"
SIMILARITIES="$3"
EDGES="$4"
NODES="$5"
THREADS="$6"
IDENTITY="$7"
BITSCORE="$8"
COVERAGE="$9"
EVALUE="${10:-1e-5}"
ANNOT_CSV="${11}"
OUT_ANNOT="${12}"

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
check_tool awk
check_tool grep
check_tool sort

# -----------------------------------------------------------------------------
# Input validation
# -----------------------------------------------------------------------------

[ -f "$INPUT" ]        || die "Input FASTA not found: $INPUT"
grep -q "^>" "$INPUT"  || die "Input does not look like a FASTA file: $INPUT"

RAW_COUNT=$(grep -c "^>" "$INPUT")
log "Input sequences: $RAW_COUNT"
[ "$RAW_COUNT" -ge 2 ] || die "Need at least 2 sequences to build a network."

# Ensure all output directories exist
for F in "$NR_FASTA" "$SIMILARITIES" "$EDGES" "$NODES"; do
    mkdir -p "$(dirname "$F")"
done

# -----------------------------------------------------------------------------
# Step 1: CD-HIT - remove redundancy
# -----------------------------------------------------------------------------

log "Step 1: CD-HIT clustering (identity=${IDENTITY})"

# CD-HIT will silently produce no output if the directory does not exist.
# Snakemake then fails because the declared output file is missing.
# The explicit mkdir -p here prevents that.
mkdir -p "$(dirname "$NR_FASTA")"

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

# -----------------------------------------------------------------------------
# Step 2: Choose search engine and run all-vs-all
# -----------------------------------------------------------------------------

BLAST_THRESHOLD=5000

if [ "$SEQ_COUNT" -lt "$BLAST_THRESHOLD" ]; then
    METHOD="blast"
else
    METHOD="mmseqs"
fi

log "Step 2: All-vs-all search using ${METHOD} (n=${SEQ_COUNT})"

# -- BLAST branch -------------------------------------------------------------
if [ "$METHOD" = "blast" ]; then

    check_tool makeblastdb
    check_tool blastp

    log "Building BLAST database..."
    makeblastdb -in "$NR_FASTA" -dbtype prot -out "${NR_FASTA}.blastdb" \
        > /dev/null 2>&1

    log "Running blastp..."
    # Columns: qseqid sseqid pident length evalue bitscore qlen slen
    #           1      2      3      4      5      6        7    8
    blastp \
        -query       "$NR_FASTA"           \
        -db          "${NR_FASTA}.blastdb" \
        -out         "$SIMILARITIES"       \
        -outfmt      "6 qseqid sseqid pident length evalue bitscore qlen slen" \
        -num_threads "$THREADS"            \
        -evalue      "$EVALUE"             \
        -seg         yes

# -- MMseqs2 branch -----------------------------------------------------------
else

    check_tool mmseqs

    TMPDIR=$(mktemp -d)
    trap 'rm -rf "$TMPDIR"' EXIT

    MMSEQS_WD="$TMPDIR/mmseqs_work"
    mkdir -p "$MMSEQS_WD"

    log "Creating MMseqs2 database..."
    mmseqs createdb "$NR_FASTA" "$TMPDIR/seqDB" > /dev/null 2>&1

    log "Running MMseqs2 search (sensitivity=7.5)..."
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
    # Columns: query target pident alnlen evalue bits qcov tcov
    #           1     2      3      4      5      6    7    8
    mmseqs convertalis \
        "$TMPDIR/seqDB"    \
        "$TMPDIR/seqDB"    \
        "$TMPDIR/resultDB" \
        "$SIMILARITIES"    \
        --format-output "query,target,pident,alnlen,evalue,bits,qcov,tcov" \
        > /dev/null 2>&1

fi

RAW_HITS=$(wc -l < "$SIMILARITIES")
log "Raw similarity hits: $RAW_HITS"
[ "$RAW_HITS" -gt 0 ] || die "No similarity hits found. Check your input or search parameters."

# -----------------------------------------------------------------------------
# Step 3: Filter edges
# -----------------------------------------------------------------------------
# BLAST:   qseqid sseqid pident alnlen evalue bitscore qlen   slen
#          1      2      3      4      5      6        7      8
#   coverage = alnlen/qlen and alnlen/slen
#
# MMseqs2: query  target pident alnlen evalue bits     qcov   tcov
#          1      2      3      4      5      6        7      8
#   coverage already provided as fractions [0-1]
# -----------------------------------------------------------------------------

log "Step 3: Filtering edges (bitscore>=${BITSCORE}, coverage>=${COVERAGE}, evalue<=${EVALUE})"

if [ "$METHOD" = "blast" ]; then
    awk \
        -v bits="$BITSCORE"   \
        -v cov="$COVERAGE"    \
        -v max_eval="$EVALUE" \
        'BEGIN { OFS="\t" }
        {
            if ($1 == $2)                 next
            if ($6 < bits)                next
            if ($5 + 0 > max_eval + 0)    next
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
            if ($1 == $2)                 next
            if ($6 < bits)                next
            if ($5 + 0 > max_eval + 0)    next
            if ($7 < cov || $8 < cov)     next
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

# -----------------------------------------------------------------------------
# Step 4: Node metadata
# Columns: id | full_header | length
# -----------------------------------------------------------------------------

log "Step 4: Extracting node metadata..."

printf "id\tfull_header\tlength\n" > "$NODES"

awk '
    /^>/ {
        if (seq_id != "") print seq_id "\t" full_header "\t" seq_len
        full_header = substr($0, 2)
        seq_id = $1; sub(/^>/, "", seq_id)
        seq_len = 0
        next
    }
    {
        gsub(/[[:space:]]/, "")
        seq_len += length($0)
    }
    END {
        if (seq_id != "") print seq_id "\t" full_header "\t" seq_len
    }
' "$NR_FASTA" >> "$NODES"

NODE_COUNT=$(( $(wc -l < "$NODES") - 1 ))
log "Nodes in metadata table: $NODE_COUNT"

# -----------------------------------------------------------------------------
# Step 5: Network statistics
# -----------------------------------------------------------------------------

log "Step 5: Network statistics"

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

# -----------------------------------------------------------------------------
# Step 6: Cytoscape-ready output
# -----------------------------------------------------------------------------

SIF="${EDGES%.tsv}.sif"
EA="${EDGES%.tsv}.ea"

log "Step 6: Writing Cytoscape SIF -> $SIF"
awk 'BEGIN{OFS="\t"} {print $1, "similarity", $2}' "$EDGES" > "$SIF"

printf "bitscore\n" > "$EA"
awk 'BEGIN{OFS="\t"} {print $1" (similarity) "$2, $3}' "$EDGES" >> "$EA"
# -----------------------------------------------------------------------------
# Step 7: Generate annotation file (class & order)
# -----------------------------------------------------------------------------

log "Step 7: Generating annotation file (class/order)..."

# Check inputs
[ -f "$ANNOT_CSV" ] || die "Annotation CSV not found: $ANNOT_CSV"

mkdir -p "$(dirname "$OUT_ANNOT")"

# Extract relevant columns dynamically
awk -F',' '
BEGIN {
    OFS="\t"
}
NR==1 {
    for(i=1;i<=NF;i++){
        if($i=="locus_tag") l=i
        if($i=="class") c=i
        if($i=="order") o=i
    }

    if(!l || !c || !o){
        print "[ERROR] Required columns (locus_tag, class, order) not found in CSV" > "/dev/stderr"
        exit 1
    }

    print "id","class","order"
    next
}
{
    id = $l
    class = $c
    order = $o

    if(id == "" || id == "NA") next

    if(class == "") class = "Unknown"
    if(order == "") order = "Unknown"

    print id, class, order
}
' "$ANNOT_CSV" \
| sort -k1,1 -u \
> "$OUT_ANNOT"

ANNOT_COUNT=$(( $(wc -l < "$OUT_ANNOT") - 1 ))
log "Annotation entries written: $ANNOT_COUNT"


# -----------------------------------------------------------------------------
# Done
# -----------------------------------------------------------------------------

log "SSN pipeline completed successfully."
log ""
log "Output files:"
log "  Non-redundant FASTA : $NR_FASTA"
log "  Raw similarities    : $SIMILARITIES"
log "  Edge list (TSV)     : $EDGES"
log "  Node metadata (TSV) : $NODES"
log "  Cytoscape SIF       : $SIF"
log "  Cytoscape EA        : $EA"
log ""
log "To visualise in Cytoscape:"
log "  1. File > Import > Network from File  -> $SIF"
log "  2. File > Import > Table from File    -> $NODES  (map 'id' to node name)"
log "  3. File > Import > Table from File    -> $EA     (edge bitscore attribute)"