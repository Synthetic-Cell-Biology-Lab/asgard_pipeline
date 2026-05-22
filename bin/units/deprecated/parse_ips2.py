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
# Step 1: Extract proteins matching the search term (any analysis)
# -------------------------------

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
# Step 2: Subset the entire interpro table to only matching proteins
# -------------------------------

print("📋 Subsetting full dataset to matching proteins...")
con.execute("""
    CREATE TABLE subset AS
    SELECT i.*
    FROM interpro i
    INNER JOIN matching_proteins m ON i.protein = m.protein
""")

n_rows = con.execute("SELECT COUNT(*) FROM subset").fetchone()[0]
print(f"  → {n_rows} rows in subset")


# -------------------------------
# Step 3: Extract Pfam domains per protein as a list
# -------------------------------

print("📐 Extracting Pfam domains per protein...")
con.execute("""
    CREATE TABLE pfam_domain_lists AS
    SELECT
        protein,
        list(sig_desc) AS domain_list
    FROM subset
    WHERE analysis = 'Pfam'
      AND sig_desc IS NOT NULL
    GROUP BY protein
    ORDER BY protein
""")

n_pfam_proteins = con.execute("SELECT COUNT(*) FROM pfam_domain_lists").fetchone()[0]
print(f"  → {n_pfam_proteins} proteins have Pfam annotations")


# -------------------------------
# Step 4: Get all unique Pfam domains for pivot columns
# -------------------------------

print("🗂️  Collecting unique Pfam domains...")
domains = [
    row[0]
    for row in con.execute("""
        SELECT DISTINCT sig_desc
        FROM subset
        WHERE analysis = 'Pfam'
          AND sig_desc IS NOT NULL
        ORDER BY sig_desc
    """).fetchall()
]

n_domains = len(domains)
print(f"  → {n_domains} unique Pfam domains found")


# -------------------------------
# Step 5: Pivot — count occurrences of each domain per protein
# (uses domain_list from Step 3 to count appearances)
# -------------------------------

print("🔄 Pivoting domain counts...")

pivot_cols = ",\n        ".join(
    f"len(list_filter(domain_list, x -> x = '{d.replace(chr(39), chr(39)*2)}')) "
    f'AS "{d}"'
    for d in domains
)

pivot_sql = f"""
    CREATE TABLE domain_matrix AS
    SELECT
        protein,
        {pivot_cols}
    FROM pfam_domain_lists
    ORDER BY protein
"""

con.execute(pivot_sql)


# -------------------------------
# Step 6: Export full matrix
# -------------------------------

print(f"📤 Exporting domain matrix → {OUTFILE}")
con.execute(f"""
    COPY (SELECT * FROM domain_matrix)
    TO '{OUTFILE}' (DELIMITER '\t', HEADER, FORMAT CSV)
""")


# -------------------------------
# Step 7: Export protein IDs only
# -------------------------------

print(f"📤 Exporting protein IDs → {PROTEIN_IDS}")
con.execute(f"""
    COPY (SELECT protein FROM domain_matrix)
    TO '{PROTEIN_IDS}' (DELIMITER '\t', HEADER FALSE, FORMAT CSV)
""")


print("✅ Done.")