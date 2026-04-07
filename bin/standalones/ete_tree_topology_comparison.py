import glob
import os
from itertools import combinations
from collections import defaultdict
import re

from Bio import SeqIO
from ete3 import Tree

TREE_DIR = "/home/anirudh/asgard_pipeline/database/protein_sets/SepF/SepF_1/final_tree_construction"
ALIGN_DIR = TREE_DIR

out_path = f"{TREE_DIR}/tree_comparison_out.txt"


def load_fasta_ids(fasta):
    return [record.id for record in SeqIO.parse(fasta, "fasta")]


def normalize_id(raw_id):
    s = raw_id.split()[0]
    s = re.sub(r'_(FtsZ1|FtsZ2|CetZ|Tubulin)$', '', s)
    s = s.replace('|', '_').replace('-', '_').replace('.', '_')
    s = re.sub(r'_+', '_', s).strip('_')
    return s


def get_basename(path):
    return os.path.basename(path).replace(".treefile", "")


with open(out_path, "w") as out:

    def log(msg):
        print(msg)
        out.write(msg + "\n")

    tree_files = glob.glob(f"{TREE_DIR}/*.treefile")
    log(f"Found {len(tree_files)} trees\n")

    for t1_file, t2_file in combinations(tree_files, 2):

        name1 = get_basename(t1_file)
        name2 = get_basename(t2_file)

        log("\n==============================")
        log(f"Comparing: {name1} vs {name2}")

        fasta1 = f"{ALIGN_DIR}/{name1}.aligned.fasta"
        fasta2 = f"{ALIGN_DIR}/{name2}.aligned.fasta"

        if not os.path.exists(fasta1) or not os.path.exists(fasta2):
            log("❌ Missing FASTA, skipping")
            continue

        # -------------------------
        # Load trees — format=100 reads IQ-TREE dual supports natively
        # -------------------------
        try:
            with open(t1_file) as f:
                t1_str = f.read().strip()
            with open(t2_file) as f:
                t2_str = f.read().strip()
            
            # Strip dual supports like 98.9/100 -> 98.9
            t1_str = re.sub(r'(\d+\.?\d*)/\d+\.?\d*', r'\1', t1_str)
            t2_str = re.sub(r'(\d+\.?\d*)/\d+\.?\d*', r'\1', t2_str)
            
            t1 = Tree(t1_str, format=1)
            t2 = Tree(t2_str, format=1)
        except Exception as e:
            log(f"❌ Tree parsing failed: {e}")
            continue

        # -------------------------
        # Find common taxa via normalized keys
        # -------------------------
        t1_key_to_leaf = defaultdict(list)
        t2_key_to_leaf = defaultdict(list)

        for leaf in t1.get_leaf_names():
            t1_key_to_leaf[normalize_id(leaf)].append(leaf)

        for leaf in t2.get_leaf_names():
            t2_key_to_leaf[normalize_id(leaf)].append(leaf)

        common_keys = set(t1_key_to_leaf.keys()) & set(t2_key_to_leaf.keys())

        log(f"Leaves t1: {len(t1.get_leaf_names())}")
        log(f"Leaves t2: {len(t2.get_leaf_names())}")
        log(f"Common taxa: {len(common_keys)}")

        if len(common_keys) < 4:
            log("⚠️ Too few common taxa, skipping")
            continue

        # -------------------------
        # Prune
        # -------------------------
        t1_keep = [t1_key_to_leaf[k][0] for k in common_keys]
        t2_keep = [t2_key_to_leaf[k][0] for k in common_keys]

        try:
            t1.prune(t1_keep)
            t2.prune(t2_keep)
        except Exception as e:
            log(f"❌ Prune failed: {e}")
            continue

        # -------------------------
        # Rename leaves to normalized keys so RF comparison works
        # -------------------------
        for node in t1.get_leaves():
            node.name = normalize_id(node.name)
        for node in t2.get_leaves():
            node.name = normalize_id(node.name)

        # -------------------------
        # RF distance
        # -------------------------
        rf, max_rf, common_leaves, parts_t1, parts_t2, *_ = t1.robinson_foulds(
            t2, unrooted_trees=True
        )

        log(f"RF distance: {rf}/{max_rf}")
        log(f"Normalized RF: {rf/max_rf:.4f}" if max_rf > 0 else "Normalized RF: N/A")

        log("\n--- Tree 1 ---")
        out.write(t1.write(format=1) + "\n\n")
        log("--- Tree 2 ---")
        out.write(t2.write(format=1) + "\n\n")