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

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

mkdir -p "$TMPDIR/cds"

# Determine the column index from the header
COL_IDX=$(head -n1 "$CSV" | tr ',' '\n' | nl -v1 | awk -v c="$COLUMN" '$2==c{print $1}')

if [[ -z "$COL_IDX" ]]; then
    echo "Column '$COLUMN' not found."
    exit 1
fi

tail -n +2 "$CSV" | while IFS=, read -r -a fields; do
    ACC=$(echo "${fields[$((COL_IDX-1))]}" | xargs)   # trim whitespace

    SRC="$DB/$ACC"

    if [[ -d "$SRC" ]]; then
        cp -a "$SRC" "$TMPDIR/"
    else
        echo "Warning: directory '$SRC' not found" >&2
    fi
done

tar --zstd -cf "${OUT}.tar.zst" -C "$TMPDIR" cds

echo "Archive written to ${OUT}.tar.zst"