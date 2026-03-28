import duckdb
import sys

PARQUET_FILE = snakemake.input.database

search_string = snakemake.params.search_string
rstring       = snakemake.params.rstring

PROTEIN_IDS = snakemake.output.protein_ids
OUTFILE     = snakemake.output.outfile


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
# Connect and load parquet
# -------------------------------

con = duckdb.connect()

print("📂 Loading parquet...")
con.execute(f"""
    CREATE TABLE interpro AS
    SELECT * FROM read_parquet('{PARQUET_FILE}')
""")


# -------------------------------
# Step 1: find proteins matching the search term (Pfam only)
# -------------------------------

print("🔍 Finding matching proteins...")
con.execute(f"""
    CREATE TABLE matching_proteins AS
    SELECT DISTINCT protein
    FROM interpro
    WHERE analysis = 'Pfam'
      AND {condition}
""")

n_proteins = con.execute("SELECT COUNT(*) FROM matching_proteins").fetchone()[0]
print(f"  → {n_proteins} proteins match the search term")

if n_proteins == 0:
    sys.exit("❌ No proteins matched the search. Check your search term.")


# -------------------------------
# Step 2: get ALL Pfam domains for those proteins
# (not just the matching domain — we want the full domain architecture)
# -------------------------------

print("📐 Fetching full Pfam domain architecture for matching proteins...")
con.execute("""
    CREATE TABLE protein_domains AS
    SELECT i.protein, i.sig_desc
    FROM interpro i
    INNER JOIN matching_proteins m ON i.protein = m.protein
    WHERE i.analysis = 'Pfam'
      AND i.sig_desc IS NOT NULL
""")


# -------------------------------
# Step 3: get the distinct domain names to pivot on
# -------------------------------

domains = [
    row[0]
    for row in con.execute(
        "SELECT DISTINCT sig_desc FROM protein_domains ORDER BY sig_desc"
    ).fetchall()
]

n_domains = len(domains)
print(f"  → {n_domains} unique Pfam domains found across matching proteins")


# -------------------------------
# Step 4: pivot — one column per domain, value = count per protein
# Each CASE counts how many times that domain appears in that protein.
# -------------------------------

pivot_cols = ",\n        ".join(
    f"COUNT(CASE WHEN sig_desc = {duckdb.typing.VARCHAR.cast(d)!r} THEN 1 END) "
    f"AS \"{d}\""
    for d in domains
)

# duckdb doesn't have a clean literal quoting helper exposed, use parameter binding
# via a formatting approach — domain names go into the SQL as string literals
pivot_cols = ",\n        ".join(
    f"COUNT(CASE WHEN sig_desc = '{d.replace(chr(39), chr(39)*2)}' THEN 1 END) "
    f'AS "{d}"'
    for d in domains
)

pivot_sql = f"""
    CREATE TABLE domain_matrix AS
    SELECT
        protein,
        {pivot_cols}
    FROM protein_domains
    GROUP BY protein
    ORDER BY protein
"""

print("🔄 Pivoting domain counts...")
con.execute(pivot_sql)


# -------------------------------
# Step 5: export full matrix
# -------------------------------

print(f"📤 Exporting domain matrix → {OUTFILE}")
con.execute(f"""
    COPY (SELECT * FROM domain_matrix)
    TO '{OUTFILE}' (DELIMITER '\t', HEADER, FORMAT CSV)
""")


# -------------------------------
# Step 6: export protein IDs only
# -------------------------------

print(f"📤 Exporting protein IDs → {PROTEIN_IDS}")
con.execute(f"""
    COPY (SELECT protein FROM domain_matrix)
    TO '{PROTEIN_IDS}' (DELIMITER '\t', HEADER FALSE, FORMAT CSV)
""")


print("✅ Done.")