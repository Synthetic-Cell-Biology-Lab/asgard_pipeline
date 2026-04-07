import os
import subprocess
import pandas as pd
from Bio import SeqIO

# =========================
# INPUTS
# =========================
CSV_FILE = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/extraction_exploration/ftsz.rev.csv"
FASTA_FILE = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/ftsz.rev.fasta"
OUT_DIR = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/mmseqs_tmp"

OUT_CSV = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/ftsz.rep.rev.csv"
OUT_FASTA = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/ftsz.rep.rev.fasta"

os.makedirs(OUT_DIR, exist_ok=True)

# =========================
# LOAD DATA
# =========================
df = pd.read_csv(CSV_FILE)
records = {rec.id: rec for rec in SeqIO.parse(FASTA_FILE, "fasta")}

# =========================
# CLASSIFY
# =========================
def classify(product):
    p = str(product).lower()
    if "ftsz1" in p:
        return "FtsZ1"
    elif "ftsz2" in p:
        return "FtsZ2"
    elif "cetz" in p:
        return "CetZ"
    elif "tubulin" in p:
        return "Tubulin"
    return None

df["family_type"] = df["Manual_annotation"].apply(classify)
df = df[df["family_type"].notna()]

# =========================
# WRITE FASTA HELPER
# =========================
def write_fasta(sub_df, path):
    recs = []
    for lid in sub_df["locus_tag"]:
        if lid in records:
            recs.append(records[lid])
    SeqIO.write(recs, path, "fasta")

# =========================
# MMSEQS MEDOID
# =========================
def mmseqs_medoid(fasta_in, tmp_dir):
    out_prefix = os.path.join(tmp_dir, "cluster")

    subprocess.run(
        f"mmseqs easy-cluster {fasta_in} {out_prefix} {tmp_dir}/tmp "
        f"--min-seq-id 0.3 -c 0.8",
        shell=True,
        check=True
    )

    rep_fasta = out_prefix + "_rep_seq.fasta"

    return [rec.id for rec in SeqIO.parse(rep_fasta, "fasta")]

# =========================
# DOWNSAMPLE PER GROUP
# =========================
def downsample_group(sub_df, label):
    selected_ids = []

    for genus, gdf in sub_df.groupby("genus"):
        genus_dir = os.path.join(OUT_DIR, f"{label}_{genus}")
        os.makedirs(genus_dir, exist_ok=True)

        fasta_path = os.path.join(genus_dir, "input.fasta")
        write_fasta(gdf, fasta_path)

        reps = mmseqs_medoid(fasta_path, genus_dir)

        selected_ids.extend(reps)

    selected_df = sub_df[sub_df["locus_tag"].isin(selected_ids)]

    # =========================
    # Ensure order coverage
    # =========================
    present_orders = set(selected_df["order"])
    all_orders = set(sub_df["order"].dropna())
    

    missing_orders = all_orders - present_orders

    for order in missing_orders:
        candidates = sub_df[sub_df["order"] == order]

        # pick one representative using mmseqs again
        tmp_dir = os.path.join(OUT_DIR, f"{label}_order_{order}")
        os.makedirs(tmp_dir, exist_ok=True)

        fasta_path = os.path.join(tmp_dir, "input.fasta")
        write_fasta(candidates, fasta_path)

        reps = mmseqs_medoid(fasta_path, tmp_dir)

        selected_ids.extend(reps)

    return sub_df[sub_df["locus_tag"].isin(selected_ids)].drop_duplicates("locus_tag")

# =========================
# MAIN LOGIC
# =========================
final_dfs = []

# Keep Tubulin
tubulin_df = df[df["family_type"] == "Tubulin"]
final_dfs.append(tubulin_df)

# Downsample others
for label in ["FtsZ1", "FtsZ2", "CetZ"]:
    sub = df[df["family_type"] == label]
    if len(sub) == 0:
        continue

    print(f"Processing {label}...")
    ds = downsample_group(sub, label)
    final_dfs.append(ds)

# Combine
final_df = pd.concat(final_dfs, ignore_index=True)

# =========================
# WRITE FASTA
# =========================
selected_ids = set(final_df["locus_tag"])
selected_records = [records[i] for i in selected_ids if i in records]

SeqIO.write(selected_records, OUT_FASTA, "fasta")
final_df.to_csv(OUT_CSV, index=False)

print("Done.")