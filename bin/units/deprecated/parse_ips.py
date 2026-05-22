import duckdb
import sys

PARQUET_FILE = snakemake.input.database

search_string = snakemake.params.search_string
rstring = snakemake.params.rstring

PROTEIN_IDS = snakemake.output.protein_ids
OUTFILE = snakemake.output.outfile


# -------------------------------
# Validate mutually exclusive options
# -------------------------------

if search_string and rstring:
    sys.exit("❌ Provide only one of 'search_string' OR 'rstring' in config.")

if not search_string and not rstring:
    sys.exit("❌ You must provide either 'search_string' OR 'rstring' in config.")


# -------------------------------
# Build SQL condition
# -------------------------------

if search_string:
    print(f"🔎 Using LIKE search: {search_string}")
    condition = f"LOWER(domains) LIKE '%{search_string.lower()}%'"

elif rstring:
    print(f"🔎 Using REGEX search: {rstring}")
    condition = f"regexp_matches(domains, '{rstring}')"


# -------------------------------
# Connect and load parquet
# -------------------------------

con = duckdb.connect()

print("📂 Loading parquet...")
con.execute(f"""
    CREATE TABLE interpro AS
    SELECT * FROM read_parquet('{PARQUET_FILE}')
""")


# -------------------------------
# Export full table
# -------------------------------

print(f"🧬 Exporting filtered domains → {OUTFILE}")

con.execute(f"""
    COPY (
        SELECT protein, domains
        FROM interpro
        WHERE domains IS NOT NULL
          AND {condition}
    )
    TO '{OUTFILE}' (DELIMITER '\t', HEADER, FORMAT CSV)
""")


# -------------------------------
# Export protein IDs only (no header)
# -------------------------------

print(f"🧬 Exporting protein IDs → {PROTEIN_IDS}")

con.execute(f"""
    COPY (
        SELECT protein
        FROM interpro
        WHERE domains IS NOT NULL
          AND {condition}
    )
    TO '{PROTEIN_IDS}' (DELIMITER '\t', HEADER FALSE, FORMAT CSV)
""")


print("✅ Done.")