import sys
print("PYTHON:", sys.executable)


import duckdb
import os
import random
import pandas as pd

SUMMARY_PARQUET = snakemake.input.database      # protein_summary parquet (ipr_desc based)
RAW_PARQUET     = snakemake.input.raw_database  # raw InterPro parquet (all analyses)

search_string = snakemake.params.search_string
rstring       = snakemake.params.rstring

OUTFILE     = snakemake.output.outfile
PROTEIN_IDS = snakemake.output.protein_ids
ITOL_DIR    = snakemake.output.itol_dir

os.makedirs(ITOL_DIR, exist_ok=True)


# -------------------------------
# Validate mutually exclusive options
# -------------------------------

if search_string and rstring:
    sys.exit("❌ Provide only one of 'search_string' OR 'rstring' in config.")

if not search_string and not rstring:
    sys.exit("❌ You must provide either 'search_string' OR 'rstring' in config.")


# -------------------------------
# Build SQL condition against the summary 'domains' column
# (same logic as the original summary-search script)
# -------------------------------

if search_string:
    print(f"🔎 Using LIKE search: {search_string}")
    condition = f"LOWER(domains) LIKE '%{search_string.lower()}%'"
elif rstring:
    print(f"🔎 Using REGEX search: {rstring}")
    condition = f"regexp_matches(domains, '{rstring}')"


# -------------------------------
# Step 1: Find matching proteins from the SUMMARY parquet
# (this is what gave you the 2500 proteins before)
# -------------------------------

con = duckdb.connect()

print("📂 Loading summary parquet...")
con.execute(f"""
    CREATE TABLE summary AS
    SELECT * FROM read_parquet('{SUMMARY_PARQUET}')
    WHERE domains IS NOT NULL
""")

print("🔍 Finding proteins matching the search term...")
con.execute(f"""
    CREATE TABLE matching_proteins AS
    SELECT DISTINCT protein
    FROM summary
    WHERE {condition}
""")

n_proteins = con.execute("SELECT COUNT(*) FROM matching_proteins").fetchone()[0]
print(f"  → {n_proteins} proteins match the search term")

if n_proteins == 0:
    sys.exit("❌ No proteins matched the search. Check your search term.")


# -------------------------------
# Step 2: Export summary TSV and protein IDs
# (same outputs as the original script)
# -------------------------------

print(f"📤 Exporting filtered summary → {OUTFILE}")
con.execute(f"""
    COPY (
        SELECT s.protein, s.domains
        FROM summary s
        INNER JOIN matching_proteins m ON s.protein = m.protein
        ORDER BY s.protein
    )
    TO '{OUTFILE}' (DELIMITER '\t', HEADER, FORMAT CSV)
""")

print(f"📤 Exporting protein IDs → {PROTEIN_IDS}")
con.execute(f"""
    COPY (SELECT protein FROM matching_proteins ORDER BY protein)
    TO '{PROTEIN_IDS}' (DELIMITER '\t', HEADER FALSE, FORMAT CSV)
""")


# -------------------------------
# Step 3: Pull ALL raw hits for matching proteins from the RAW parquet
# This is the key fix — we use the summary to find proteins,
# then go back to the raw data to get every database hit for iTOL
# -------------------------------

print("📂 Loading raw InterPro parquet...")
con.execute(f"""
    CREATE TABLE raw AS
    SELECT * FROM read_parquet('{RAW_PARQUET}')
""")

print("📋 Subsetting raw hits to matching proteins (all databases)...")
con.execute("""
    CREATE TABLE subset AS
    SELECT r.*
    FROM raw r
    INNER JOIN matching_proteins m ON r.protein = m.protein
""")

n_rows = con.execute("SELECT COUNT(*) FROM subset").fetchone()[0]
print(f"  → {n_rows} raw rows across all databases")

df = con.execute("SELECT * FROM subset").df()

databases = sorted(df["analysis"].unique())
print(f"  → {len(databases)} databases found: {', '.join(databases)}")


# -------------------------------
# Aesthetic color map — vivid, well-separated colors per sig_acc
# Fixed seed for reproducibility
# -------------------------------

rng = random.Random(42)

def rand_vivid_color(rng):
    """One channel high, one low, one mid — always a vivid hue."""
    channels = [rng.randint(180, 255), rng.randint(0, 75), rng.randint(75, 180)]
    rng.shuffle(channels)
    return "#{:02x}{:02x}{:02x}".format(*channels)

all_accs  = df[df["sig_acc"].notna()]["sig_acc"].unique()
color_map = {acc: rand_vivid_color(rng) for acc in all_accs}


# -------------------------------
# Shape map — one distinct iTOL shape per database
# -------------------------------

SHAPES   = ["EL", "HH", "OC", "DI", "TR", "TL", "PU", "RE"]
db_shape = {db: SHAPES[i % len(SHAPES)] for i, db in enumerate(databases)}

print("  → Shape assignments:")
for db, shape in db_shape.items():
    print(f"      {db}: {shape}")


# -------------------------------
# Write one iTOL DATASET_DOMAINS file per database
# -------------------------------

def write_itol(df_full, db_name, shape, out_dir, color_map):
    all_proteins = sorted(df_full["protein"].unique())

    df_db = df_full[
        (df_full["analysis"] == db_name) &
        (df_full["sig_acc"].notna())
    ]
    db_proteins = set(df_db["protein"].unique())

    lines   = []
    skipped = 0

    for protein in all_proteins:
        protein_rows = df_full[df_full["protein"] == protein]
        length_vals  = protein_rows["length"].dropna()

        if length_vals.empty:
            skipped += 1
            continue

        length = int(length_vals.iloc[0])

        if protein not in db_proteins:
            # No hits in this database — emit backbone-only so protein stays visible
            lines.append(f"{protein},{length}")
            continue

        group   = df_db[df_db["protein"] == protein].sort_values("start")
        domains = []
        for _, row in group.iterrows():
            label = row["sig_desc"] if pd.notna(row["sig_desc"]) else row["sig_acc"]
            label = str(label).replace("|", " ").replace(",", ";")
            color = color_map[row["sig_acc"]]
            domains.append(f"{shape}|{int(row['start'])}|{int(row['end'])}|{color}|{label}")

        lines.append(f"{protein},{length}," + ",".join(domains))

    outfile = os.path.join(out_dir, f"itol_{db_name.lower()}.txt")
    with open(outfile, "w") as f:
        f.write("DATASET_DOMAINS\n")
        f.write("SEPARATOR COMMA\n")
        f.write(f"DATASET_LABEL,{db_name}_domains\n")
        f.write("COLOR,#ff0000\n")
        f.write("SHOW_DOMAIN_LABELS,1\n")
        f.write("LABEL_SIZE_FACTOR,0.8\n")
        f.write("LABEL_AUTO_COLOR,1\n")
        f.write("BACKBONE_COLOR,#cccccc\n")
        f.write("BACKBONE_HEIGHT,6\n\n")
        f.write("DATA\n")
        f.write("\n".join(lines))

    with_domains = sum(1 for l in lines if "|" in l)
    print(f"  [DONE] {outfile}  ({len(lines)} proteins, {with_domains} with {db_name} domains)")
    if skipped:
        print(f"    ⚠️  {skipped} proteins skipped — no length info anywhere in raw data")


print(f"🎨 Writing iTOL domain files → {ITOL_DIR}")
for db in databases:
    write_itol(df, db, db_shape[db], ITOL_DIR, color_map)


print("✅ Done.")