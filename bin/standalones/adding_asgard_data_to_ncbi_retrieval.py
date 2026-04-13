import pandas as pd

# Load files
df1 = pd.read_csv("/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/extraction_exploration/ftsz_csm_clean.rev.csv")  # seq_id file (your FtsZ/tubulin list)
df2 = pd.read_csv("/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/extraction_exploration/ftsz.rev.csv")  # annotation file

# ---- Clean column names ----
df1.columns = df1.columns.str.strip()
df2.columns = df2.columns.str.strip()

# ---- Ensure matching keys ----
df1["seq_id"] = df1["seq_id"].astype(str).str.strip()
df2["locus_tag"] = df2["locus_tag"].astype(str).str.strip()

# ---- Select only useful columns from file2 ----
cols_to_add = [
    "locus_tag",
    "product",
    "gene",
    "protein_length",
    "Organism Name",
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
]

df2_sub = df2[cols_to_add]

# ---- Merge ----
merged = pd.merge(
    df1,
    df2_sub,
    left_on="seq_id",
    right_on="locus_tag",
    how="left"   # keep all from file1
)

# ---- Drop redundant column ----
merged = merged.drop(columns=["locus_tag"])

# ---- Save ----
merged.to_csv("/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/extraction_exploration/ftsz_csm_clean_fin.rev.csv", index=False)

print("Done. Rows:", len(merged))