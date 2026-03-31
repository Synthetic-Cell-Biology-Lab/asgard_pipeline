#!/usr/bin/env bash

set -euo pipefail

INPUT_FASTA=$1
PREFIX=$2
THREADS=${3:-8}

if [ ! -f "$INPUT_FASTA" ]; then
    echo "❌ Input FASTA not found: $INPUT_FASTA"
    exit 1
fi

echo "=============================================="
echo "🌳 Exploratory FastTree"
echo "Input: $INPUT_FASTA"
echo "Prefix: $PREFIX"
echo "Threads: $THREADS"
echo "=============================================="

############################################
# 1️⃣ Alignment (MAFFT)
############################################

echo "🧬 Running MAFFT..."

echo "Running FAMSA2 alignment..."

famsa -t "$THREADS" \
      "$INPUT_FASTA" \
      "${PREFIX}.aligned.fasta" 

echo "Alignment complete → ${PREFIX}.aligned.fasta"

############################################
# 2️⃣ FastTree
############################################

echo "🌲 Running FastTree..."

export OMP_NUM_THREADS=$THREADS

FastTreeMP \
    -lg \
    -gamma \
    "${PREFIX}.aligned.fasta" \
    > "${PREFIX}.treefile"

echo "✅ Tree written → ${PREFIX}.treefile"
echo "=============================================="