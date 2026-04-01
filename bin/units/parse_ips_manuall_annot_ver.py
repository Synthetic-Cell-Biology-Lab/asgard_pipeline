#!/usr/bin/env python3

import duckdb
import sys
import os
import random
import pandas as pd

# ============================================================
# INPUTS
# ============================================================

RAW_PARQUET = snakemake.input.raw_database
IDS_FILE    = snakemake.input.protein_ids

OUTFILE  = snakemake.output.outfile   # TSV summary
ITOL_DIR = snakemake.output.itol_dir

os.makedirs(ITOL_DIR, exist_ok=True)

# ============================================================
# LOAD IDS
# ============================================================

print("📂 Loading protein ID list...")
with open(IDS_FILE) as f:
    protein_ids = set(
        line.strip().replace(" ", "_")
        for line in f if line.strip()
    )

print(f"  → {len(protein_ids)} protein IDs loaded")

if not protein_ids:
    sys.exit("❌ No protein IDs provided")

# ============================================================
# CONNECT DB
# ============================================================

con = duckdb.connect()

# ============================================================
# CREATE ID TABLE
# ============================================================

con.execute("CREATE TABLE id_list(protein VARCHAR)")
con.executemany(
    "INSERT INTO id_list VALUES (?)",
    [(pid,) for pid in protein_ids]
)

# ============================================================
# LOAD RAW + FILTER
# ============================================================

print("📂 Loading raw InterPro parquet...")
con.execute(f"""
    CREATE TABLE raw AS
    SELECT * FROM read_parquet('{RAW_PARQUET}')
""")

print("🔍 Filtering raw data using ID list...")
con.execute("""
    CREATE TABLE subset AS
    SELECT r.*
    FROM raw r
    INNER JOIN id_list i ON r.protein = i.protein
""")

n_rows = con.execute("SELECT COUNT(*) FROM subset").fetchone()[0]
print(f"  → {n_rows} rows retained")

if n_rows == 0:
    sys.exit("❌ No matching proteins found in raw data")

# ============================================================
# EXPORT TSV (replacement for summary)
# ============================================================

print(f"📤 Writing TSV summary → {OUTFILE}")

df_subset = con.execute("SELECT * FROM subset").df()

# collapse domains per protein
df_tsv = (
    df_subset[df_subset["sig_acc"].notna()]
    .groupby("protein")
    .apply(lambda x: "; ".join(
        sorted(set(
            str(v) for v in x["sig_desc"].dropna()
        ))
    ))
    .reset_index(name="domains")
)

df_tsv.to_csv(OUTFILE, sep="\t", index=False)

# ============================================================
# PREP FOR ITOL
# ============================================================

df = df_subset
databases = sorted(df["analysis"].dropna().unique())

print(f"  → {len(databases)} databases: {', '.join(databases)}")

# ============================================================
# COLOR MAP
# ============================================================

rng = random.Random(42)

def rand_vivid_color(rng):
    channels = [rng.randint(180, 255), rng.randint(0, 75), rng.randint(75, 180)]
    rng.shuffle(channels)
    return "#{:02x}{:02x}{:02x}".format(*channels)

all_accs  = df[df["sig_acc"].notna()]["sig_acc"].unique()
color_map = {acc: rand_vivid_color(rng) for acc in all_accs}

# ============================================================
# SHAPES
# ============================================================

SHAPES   = ["EL", "HH", "OC", "DI", "TR", "TL", "PU", "RE"]
db_shape = {db: SHAPES[i % len(SHAPES)] for i, db in enumerate(databases)}

# ============================================================
# ITOL WRITER
# ============================================================

def write_itol(df_full, db_name, shape, out_dir, color_map):

    all_proteins = sorted(df_full["protein"].unique())

    df_db = df_full[
        (df_full["analysis"] == db_name) &
        (df_full["sig_acc"].notna())
    ]

    db_proteins = set(df_db["protein"].unique())

    lines = []

    for protein in all_proteins:

        protein_rows = df_full[df_full["protein"] == protein]
        length_vals  = protein_rows["length"].dropna()

        if length_vals.empty:
            continue

        length = int(length_vals.iloc[0])

        if protein not in db_proteins:
            lines.append(f"{protein},{length}")
            continue

        group = df_db[df_db["protein"] == protein].sort_values("start")

        domains = []
        for _, row in group.iterrows():
            label = row["sig_desc"] if pd.notna(row["sig_desc"]) else row["sig_acc"]
            label = str(label).replace("|", " ").replace(",", ";")

            color = color_map[row["sig_acc"]]

            domains.append(
                f"{shape}|{int(row['start'])}|{int(row['end'])}|{color}|{label}"
            )

        lines.append(f"{protein},{length}," + ",".join(domains))

    outfile = os.path.join(out_dir, f"itol_{db_name.lower()}.txt")

    with open(outfile, "w") as f:
        f.write("DATASET_DOMAINS\n")
        f.write("SEPARATOR COMMA\n")
        f.write(f"DATASET_LABEL,{db_name}_domains\n")
        f.write("COLOR,#ff0000\n")
        f.write("SHOW_DOMAIN_LABELS,1\n")
        f.write("BACKBONE_COLOR,#cccccc\n")
        f.write("BACKBONE_HEIGHT,6\n\n")
        f.write("DATA\n")
        f.write("\n".join(lines))

    print(f"  [DONE] {outfile}")

# ============================================================
# GENERATE ITOL FILES
# ============================================================

print(f"🎨 Writing iTOL files → {ITOL_DIR}")

for db in databases:
    write_itol(df, db, db_shape[db], ITOL_DIR, color_map)

print("✅ Done.")