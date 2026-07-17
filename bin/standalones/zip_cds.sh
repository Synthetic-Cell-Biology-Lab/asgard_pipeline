#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <csv> <accession_column> <cds_database> <output_prefix>"
    exit 1
fi

CSV="$1"
COLUMN="$2"
DB="$3"
OUT="$4"

# Extract accessions with a real CSV parser (handles quoting/commas correctly)
mapfile -t ACCESSIONS < <(python3 - "$CSV" "$COLUMN" <<'EOF'
import csv, sys
path, column = sys.argv[1], sys.argv[2]
with open(path, newline='') as f:
    reader = csv.DictReader(f)
    if column not in reader.fieldnames:
        sys.exit(f"Column '{column}' not found. Available: {reader.fieldnames}")
    for row in reader:
        acc = (row[column] or "").strip()
        if acc:
            print(acc)
EOF
)

if [[ ${#ACCESSIONS[@]} -eq 0 ]]; then
    echo "Error: no accessions extracted from '$CSV' column '$COLUMN'." >&2
    exit 1
fi

copied=0
missing=0
FILES_LIST=$(mktemp)
trap 'rm -f "$FILES_LIST"' EXIT


for ACC in "${ACCESSIONS[@]}"; do
    if [[ -d "$DB/$ACC" ]]; then
        echo "$ACC" >> "$FILES_LIST"
        copied=$((copied+1))
    else
        echo "Warning: directory '$DB/$ACC' not found" >&2
        missing=$((missing+1))
    fi
done


echo "Found: $copied   Missing: $missing"

if [[ "$copied" -eq 0 ]]; then
    echo "Error: no genome directories found — check DB path and column name." >&2
    exit 1
fi


tar --zstd \
    -C "$DB" \
    --transform 's,^,cds/,' \
    -cf "${OUT}.tar.zst" \
    --files-from="$FILES_LIST"

echo "Archive written to ${OUT}.tar.zst"