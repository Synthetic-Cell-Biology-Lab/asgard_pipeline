#!/usr/bin/env bash
set -euo pipefail

INPUT=$1
NR_FASTA=$2
SIMILARITIES=$3
EDGES=$4
NODES=$5
THREADS=$6
IDENTITY=$7
BITSCORE=$8
COVERAGE=$9

echo "Step 1: CD-HIT"

cd-hit \
-i "$INPUT" \
-o "$NR_FASTA" \
-c "$IDENTITY" \
-n 5 \
-T "$THREADS"

echo "Step 2: Count sequences"

SEQ_COUNT=$(grep -c "^>" "$NR_FASTA")
echo "Sequences after CD-HIT: $SEQ_COUNT"

if [ "$SEQ_COUNT" -lt 5000 ]; then
    METHOD="blast"
else
    METHOD="mmseqs"
fi

echo "Using search engine: $METHOD"

if [ "$METHOD" = "blast" ]; then

    makeblastdb -in "$NR_FASTA" -dbtype prot

    blastp \
    -query "$NR_FASTA" \
    -db "$NR_FASTA" \
    -out "$SIMILARITIES" \
    -outfmt "6 qseqid sseqid pident length evalue bitscore qlen slen"

else

    TMPDIR=$(mktemp -d)

    mmseqs createdb "$NR_FASTA" seqDB

    mmseqs search seqDB seqDB resultDB "$TMPDIR" \
    --threads "$THREADS"

    mmseqs convertalis seqDB seqDB resultDB \
    "$SIMILARITIES" \
    --format-output "query,target,pident,alnlen,evalue,bits,qcov,tcov"

    rm -rf "$TMPDIR"

fi

echo "Step 3: Filtering edges"

awk -v bits="$BITSCORE" -v cov="$COVERAGE" '{
    if ($1 != $2 && $6 >= bits && $7 >= cov && $8 >= cov)
        print $1"\t"$2"\t"$6
}' "$SIMILARITIES" > "$EDGES"

echo "Step 4: Node metadata"

grep "^>" "$NR_FASTA" \
| sed 's/>//' \
> "$NODES"

echo "SSN pipeline completed"