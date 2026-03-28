


set -euo pipefail

# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <nr.fasta> <nodes.tsv>"
    exit 1
fi




NR_FASTA="$1"
NODES="$2"




log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# -----------------------------------------------------------------------------
# Node metadata
# Columns: id | full_header | length
# Parsed from the NR FASTA — one row per sequence
# -----------------------------------------------------------------------------




log "Extracting node metadata from $NR_FASTA..."

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
