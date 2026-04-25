#!/usr/bin/env bash

set -euo pipefail

############################################
# Usage
############################################
# bash phylogeny_pipeline.sh input.fasta output_prefix threads
#
# Example:
# bash phylogeny_pipeline.sh ftsz.rev.fasta ftsz_phylo 16
############################################

INPUT_FASTA=$1
PREFIX=$2
THREADS=${3:-8}

############################################
# Check inputs
############################################

if [ ! -f "$INPUT_FASTA" ]; then
    echo "❌ Input FASTA not found: $INPUT_FASTA"
    exit 1
fi

echo "=============================================="
echo "PHYLOGENY PIPELINE STARTED"
echo "Input: $INPUT_FASTA"
echo "Prefix: $PREFIX"
echo "Threads: $THREADS"
echo "Start Time: $(date)"
echo "=============================================="

# ############################################
# # 1️⃣ Alignment (MAFFT)
# ############################################

echo "🧬 Running MAFFT alignment..."

mafft --localpair --maxiterate 1000 --thread "$THREADS" "$INPUT_FASTA" \
    > "${PREFIX}.aligned.fasta"

echo "✅ Alignment complete → ${PREFIX}.aligned.fasta"

############################################
# 1️⃣ Alignment (FAMSA2)
############################################

# echo "Running FAMSA2 alignment..."

# famsa -t "$THREADS" \
#       "$INPUT_FASTA" \
#       "${PREFIX}.aligned.fasta" 

# echo "Alignment complete → ${PREFIX}.aligned.fasta"


############################################
# 2️⃣ Trimming (ClipKIT)
############################################
# Preferred over TrimAl for protein datasets

echo "Running ClipKIT trimming..."

clipkit "${PREFIX}.aligned.fasta" \
    -m smart-gap \
    -o "${PREFIX}.trimmed.fasta"

echo "Trimming complete → ${PREFIX}.trimmed.fasta"

############################################
# 3️⃣ Model Selection + 4️⃣ ML Tree (IQ-TREE3)
############################################

echo "Running IQ-TREE3 (Model selection + ML tree)..."

# Step 1 — tree inference
iqtree3 \
    -s "${PREFIX}.trimmed.fasta" \
    -T AUTO \
    -m LG+C40+F+R \
    -bb 1000 -bnni -nstop 200 \
    -alrt 1000 \
    --runs 5 \
    -redo \
    -pre "${PREFIX}" \
    2>&1 | tee "${PREFIX}.iqtree_console.log"

# Step 2 — sCF annotation on the resulting tree
iqtree3 \
    -t "${PREFIX}.treefile" \
    -s "${PREFIX}.trimmed.fasta" \
    --scf 100 \
    -pre "${PREFIX}.cf" \
    -T AUTO \
    2>&1 | tee "${PREFIX}.scf.log"
    
echo "✅ Phylogeny complete"

############################################
# Summary
############################################


echo "=============================================="
echo "PHYLOGENY PIPELINE COMPLETE"
echo "Final Tree: ${PREFIX}.treefile"
echo "Bootstrap Tree: ${PREFIX}.contree"
echo "Log File: ${PREFIX}.log"
echo "Model Info: ${PREFIX}.iqtree"
echo "End Time: $(date)"
echo "=============================================="