set -euo pipefail

TIMESTAMP=$(date +"%Y-%m-%d_%H.%M")

QUERY=${1:-}
VERSION="1"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

QUERY_FILE="${BASE_DIR}/database/Blast/queries/${QUERY}.fasta"

OUTPUT_DIR="${BASE_DIR}/database/Blast/results/${QUERY}"
mkdir -p $OUTPUT_DIR

OUTPUT_FILE="$OUTPUT_DIR/${QUERY}_Version${VERSION}_DB_${TIMESTAMP}.tsv"
LOGFILE="$OUTPUT_DIR/${QUERY}_Version${VERSION}_DB_${TIMESTAMP}.log"


if [ -z "$QUERY_FILE" ]; then
  echo "Usage: ./run_blast.sh <query.fasta>"
  exit 1
fi

if [ ! -f "$QUERY_FILE" ]; then
  echo "❌ Query file not found: $QUERY"
  exit 1
fi


TARGET_DB="${BASE_DIR}/database/collated/Version${VERSION}/filtered/85comp10con/fasta/v${VERSION}_cp85_con10.fasta"

log() {
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] $1" | tee -a "$LOGFILE"
}

log "Starting BLASTp Search for ${QUERY}"
log "Target DB       : ${TARGET_DB}"
log "Output Directory: ${OUTPUT_DIR}"

log "Checking input file formats"
file "$TARGET_DB" >> "$LOGFILE" 2>&1
file "$QUERY_FILE" >> "$LOGFILE" 2>&1


log "Building BLAST protein database"
makeblastdb -in "$TARGET_DB" -dbtype prot >> "$LOGFILE" 2>&1
log "BLAST database construction completed"


# -----------------------------
# Run BLASTp
# -----------------------------
log "Running BLASTp search (e-value = 1e-10)"

HEADER="qseqid\tsseqid\tpident\tlength\tqlen\tslen\tevalue\tbitscore"

{
  echo -e "$HEADER"
  blastp \
    -query "$QUERY_FILE" \
    -db "$TARGET_DB" \
    -evalue 1e-10 \
    -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore"
} > "$OUTPUT_FILE" 2>> "$LOGFILE"

log "BLASTp search completed"

# -----------------------------
# Extract candidate sequences
# -----------------------------
log "Extracting full-length actin candidate protein sequences"

awk 'NR>1 {print $2}' "$OUTPUT_FILE" | sort -u \
  > "$OUTPUT_DIR/${QUERY}_candidate_ids_${TIMESTAMP}.txt"

seqkit grep \
  -f "$OUTPUT_DIR/${QUERY}_candidate_ids_${TIMESTAMP}.txt" \
  "$TARGET_DB" \
  > "$OUTPUT_DIR/${QUERY}_candidates_${TIMESTAMP}.faa"

log "Protein sequence extraction completed"

# -----------------------------
# End
# -----------------------------
log "Pipeline finished successfully"
log "Results saved in: $OUTPUT_DIR"
