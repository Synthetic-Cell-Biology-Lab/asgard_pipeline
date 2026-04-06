#!/bin/bash

# =========================
# CONFIG
# =========================
TREE_DIR="/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/with_csm_seq"              # folder with .nex files
# ALIGNMENT="alignment.fasta"  # required for AU test
OUTDIR="/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/with_csm_seq/tree_comparisons"

mkdir -p "$OUTDIR/nwk"
mkdir -p "$OUTDIR/rf"
# mkdir -p "$OUTDIR/au"

# =========================
# STEP 1: Collect treefiles
# =========================
echo "Collecting treefiles..."

trees=("$TREE_DIR/"*.treefile)
n=${#trees[@]}

echo "Found $n trees"

echo "Running RF comparisons..."

for ((i=0; i<n; i++)); do
    for ((j=i+1; j<n; j++)); do

        t1=${trees[$i]}
        t2=${trees[$j]}

        name1=$(basename "$t1" .treefile)
        name2=$(basename "$t2" .treefile)

        echo "RF: $name1 vs $name2"

        iqtree3 -rf "$t1" "$t2" > "$OUTDIR/rf/${name1}_vs_${name2}.txt"
    done
done


# # =========================
# # STEP 4: All-vs-all AU tests
# # =========================

# echo "Running AU tests..."

# for ((i=0; i<n; i++)); do
#     for ((j=i+1; j<n; j++)); do

#         t1=${trees[$i]}
#         t2=${trees[$j]}

#         name1=$(basename "$t1" .nwk)
#         name2=$(basename "$t2" .nwk)

#         pairfile="$OUTDIR/au/${name1}_vs_${name2}.trees"

#         # Combine trees
#         cat "$t1" "$t2" > "$pairfile"

#         echo "AU: $name1 vs $name2"

#         iqtree2 -s "$ALIGNMENT" \
#                  -z "$pairfile" \
#                  -zb 10000 \
#                  -au \
#                  -pre "$OUTDIR/au/${name1}_vs_${name2}"
#     done
# done

# echo "Done."