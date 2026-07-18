# python script to take a parquet file and filter by the acc column of the parquet
# using pfam ids/database ids given in config file

import pandas as pd
import duckdb
import argparse
import random
import os

# -------------------------------
# Write one combined iTOL DATASET_DOMAINS file (all databases together)
# -------------------------------

# -------------------------------
# Write one iTOL DATASET_DOMAINS file per database
# -------------------------------


def write_itol(df_full, db_name, shape, out_dir, color_map):
    all_proteins = sorted(df_full["protein"].unique())

    df_db = df_full[(df_full["analysis"] == db_name) & (df_full["sig_acc"].notna())]
    db_proteins = set(df_db["protein"].unique())

    lines = []
    skipped = 0

    for protein in all_proteins:
        protein_rows = df_full[df_full["protein"] == protein]
        length_vals = protein_rows["length"].dropna()

        if length_vals.empty:
            skipped += 1
            continue

        length = int(length_vals.iloc[0])

        if protein not in db_proteins:
            # No hits in this database — emit backbone-only so protein stays visible
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
        f.write("LABEL_SIZE_FACTOR,0.8\n")
        f.write("LABEL_AUTO_COLOR,1\n")
        f.write("BACKBONE_COLOR,#cccccc\n")
        f.write("BACKBONE_HEIGHT,6\n\n")
        f.write("DATA\n")
        f.write("\n".join(lines))

    with_domains = sum(1 for l in lines if "|" in l)
    print(
        f"  [DONE] {outfile}  ({len(lines)} proteins, {with_domains} with {db_name} domains)"
    )
    if skipped:
        print(
            f"    ⚠️  {skipped} proteins skipped — no length info anywhere in raw data"
        )


def main():
    parser = argparse.ArgumentParser()

    # --- inputs ---
    parser.add_argument(
        "--acc", nargs="+", required=True, help="List of 1 or more accession ids"
    )
    parser.add_argument(
        "-p",
        "--parquet",
        required=True,
        help="The parquet file containing interproscan results",
    )

    # --- outputs (map 1:1 to `output:` in the Snakemake rule) ---
    parser.add_argument(
        "--protein-ids",
        required=True,
        help="Output path: protein_ids.txt",
    )
    parser.add_argument(
        "--matching-tsv",
        required=True,
        help="Output path: matching_hits.tsv",
    )
    parser.add_argument(
        "--itol",
        required=True,
        help="Output path: Directory path for output ITOL domain files",
    )

    args = parser.parse_args()

    con = duckdb.connect()

    # Snakemake normally creates output dirs itself, but this keeps the
    # script safe to run standalone too.
    for path in (args.protein_ids, args.matching_tsv, args.itol):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    os.makedirs(args.itol, exist_ok=True)

    print("📂 Loading raw InterPro parquet...")
    con.execute(f"""
        CREATE TABLE raw AS
        SELECT * FROM read_parquet('{args.parquet}')
    """)
    pfam = [acc for acc in args.acc if acc.startswith("PF")]
    ipr = [acc for acc in args.acc if acc.startswith("IPR")]
    clauses = []
    params = []

    if pfam:
        clauses.append(f"sig_acc IN ({','.join(['?'] * len(pfam))})")
        params.extend(pfam)

    if ipr:
        clauses.append(f"ipr_acc IN ({','.join(['?'] * len(ipr))})")
        params.extend(ipr)

    where_clauses = " OR ".join(clauses)

    print("🔍 Finding proteins matching the search term(s)...")
    con.execute(
        f"""
        CREATE TABLE matching_proteins AS
        SELECT DISTINCT protein
        FROM raw
        WHERE {where_clauses}
        """,
        params,
    )

    con.execute(f"""
        COPY (
            SELECT protein
            FROM matching_proteins
            ORDER BY protein
        )
        TO '{args.protein_ids}'
        (FORMAT CSV, DELIMITER '\t', HEADER FALSE);
    """)
    print(f"✓ Wrote protein list -> {args.protein_ids}")

    print("📋 Subsetting raw hits to matching proteins (all databases)...")
    con.execute("""
        CREATE TABLE subset AS
        SELECT r.*
        FROM raw r
        INNER JOIN matching_proteins m ON r.protein = m.protein
    """)

    n_rows = con.execute("SELECT COUNT(*) FROM subset").fetchone()[0]
    print(f"  → {n_rows} raw rows across all databases")

    con.execute(f"""
        COPY (
            SELECT *
            FROM subset
        )
        TO '{args.matching_tsv}'
        (FORMAT CSV, DELIMITER '\t', HEADER TRUE);
    """)
    print(f"✓ Wrote matching hits -> {args.matching_tsv}")

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

    all_accs = df[df["sig_acc"].notna()]["sig_acc"].unique()
    color_map = {acc: rand_vivid_color(rng) for acc in all_accs}

    # -------------------------------
    # Shape map — one distinct iTOL shape per database
    # -------------------------------
    SHAPES = ["EL", "HH", "OC", "DI", "TR", "TL", "PU", "RE"]
    db_shape = {db: SHAPES[i % len(SHAPES)] for i, db in enumerate(databases)}

    print("  → Shape assignments:")
    for db, shape in db_shape.items():
        print(f"      {db}: {shape}")

    print("\n🖍 Writing combined iTOL dataset...")
    for db in databases:
        write_itol(df, db, db_shape[db], args.itol, color_map)


if __name__ == "__main__":
    main()
