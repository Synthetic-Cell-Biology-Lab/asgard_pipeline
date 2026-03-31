import duckdb
import sys
import os
import random
import pandas as pd

PARQUET_FILE = snakemake.input.database

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
# Build SQL condition on sig_desc
# -------------------------------

if search_string:
    print(f"🔎 Using LIKE search: {search_string}")
    condition = f"LOWER(sig_desc) LIKE '%{search_string.lower()}%'"
elif rstring:
    print(f"🔎 Using REGEX search: {rstring}")
    condition = f"regexp_matches(sig_desc, '{rstring}')"


# -------------------------------
# Connect, load parquet, filter to matching proteins
# -------------------------------

con = duckdb.connect()

print("📂 Loading parquet...")
con.execute(f"""
    CREATE TABLE interpro AS
    SELECT * FROM read_parquet('{PARQUET_FILE}')
""")

print("🔍 Finding proteins matching the search term...")
con.execute(f"""
    CREATE TABLE matching_proteins AS
    SELECT DISTINCT protein
    FROM interpro
    WHERE {condition}
""")

n_proteins = con.execute("SELECT COUNT(*) FROM matching_proteins").fetchone()[0]
print(f"  → {n_proteins} proteins match the search term")

if n_proteins == 0:
    sys.exit("❌ No proteins matched the search. Check your search term.")


# -------------------------------
# Subset full table to matching proteins (all databases)
# -------------------------------

print("📋 Subsetting full dataset to matching proteins (all databases)...")
con.execute("""
    CREATE TABLE subset AS
    SELECT i.*
    FROM interpro i
    INNER JOIN matching_proteins m ON i.protein = m.protein
""")

n_rows = con.execute("SELECT COUNT(*) FROM subset").fetchone()[0]
print(f"  → {n_rows} rows in subset")


# -------------------------------
# Export full subset as TSV
# -------------------------------

print(f"📤 Exporting all hits → {OUTFILE}")
con.execute(f"""
    COPY (SELECT * FROM subset ORDER BY protein)
    TO '{OUTFILE}' (DELIMITER '\t', HEADER, FORMAT CSV)
""")


# -------------------------------
# Export protein IDs only
# -------------------------------

print(f"📤 Exporting protein IDs → {PROTEIN_IDS}")
con.execute(f"""
    COPY (SELECT DISTINCT protein FROM subset ORDER BY protein)
    TO '{PROTEIN_IDS}' (DELIMITER '\t', HEADER FALSE, FORMAT CSV)
""")


# -------------------------------
# Load subset into pandas for iTOL generation
# -------------------------------

print("🐼 Loading subset into pandas...")
df = con.execute("SELECT * FROM subset").df()

# NOTE: do NOT filter out null sig_acc at the dataframe level —
# that would drop entire proteins that only have hits in databases
# like Gene3D or Coils where sig_acc can be null for some rows.
# Instead, null sig_acc rows are skipped at the per-domain level below.

databases = sorted(df["analysis"].unique())
print(f"  → {len(databases)} databases found: {', '.join(databases)}")


# -------------------------------
# Aesthetic color map — vivid, well-separated colors per sig_acc
# Fixed seed for reproducibility across runs
# -------------------------------

rng = random.Random(42)

def rand_vivid_color(rng):
    """One channel high (180-255), one low (0-75), one mid (75-180)
    — guarantees a vivid hue rather than muddy grey."""
    channels = [rng.randint(180, 255), rng.randint(0, 75), rng.randint(75, 180)]
    rng.shuffle(channels)
    return "#{:02x}{:02x}{:02x}".format(*channels)

all_accs  = df[df["sig_acc"].notna()]["sig_acc"].unique()
color_map = {acc: rand_vivid_color(rng) for acc in all_accs}


# -------------------------------
# Shape map — one distinct iTOL shape per database
# RE=rectangle, EL=ellipse, HH=horiz hexagon, OC=octagon,
# DI=diamond, TR=right triangle, TL=left triangle, PU=up pentagram
# -------------------------------

SHAPES = ["EL", "HH", "OC", "DI", "TR", "TL", "PU", "RE"]

db_shape = {db: SHAPES[i % len(SHAPES)] for i, db in enumerate(databases)}
print("  → Shape assignments:")
for db, shape in db_shape.items():
    print(f"      {db}: {shape}")


# -------------------------------
# Write one iTOL DATASET_DOMAINS file per database
# -------------------------------

def write_itol(df_full, db_name, shape, out_dir, color_map):
    """
    df_full  : the complete subset dataframe (all databases)
    db_name  : the analysis name for this file
    shape    : iTOL shape string for this database
    """
    # All proteins that appear anywhere in the subset — ensures no protein
    # is dropped just because it lacks hits in this specific database.
    all_proteins = sorted(df_full["protein"].unique())

    # Per-database hits (only rows with a valid sig_acc are drawable)
    df_db = df_full[
        (df_full["analysis"] == db_name) &
        (df_full["sig_acc"].notna())
    ]
    db_proteins = set(df_db["protein"].unique())

    lines = []
    skipped = 0

    for protein in all_proteins:
        # Get length from any row for this protein across all databases
        protein_rows = df_full[df_full["protein"] == protein]
        length_vals  = protein_rows["length"].dropna()

        if length_vals.empty:
            skipped += 1
            continue

        length = int(length_vals.iloc[0])

        if protein not in db_proteins:
            # Protein matched the search but has no hits in this database —
            # still emit a backbone-only line so it appears in the tree
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
        print(f"    ⚠️  {skipped} proteins skipped — no length info in any database row")


print(f"🎨 Writing iTOL domain files → {ITOL_DIR}")
for db in databases:
    write_itol(df, db, db_shape[db], ITOL_DIR, color_map)


print("✅ Done.")