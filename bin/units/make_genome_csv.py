import traceback

import pandas as pd


# ---------------------------------------------------------
# Snakemake inputs / outputs / params / log
# ---------------------------------------------------------

checkm2_file = snakemake.input.checkm2
taxonomy_file = snakemake.input.taxonomy

output_unfiltered = snakemake.output.unfiltered
output_filtered = snakemake.output.filtered

min_completeness = float(snakemake.params.completeness)
max_contamination = float(snakemake.params.contamination)

log_path = snakemake.log[0]

TAX_RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]

CHECKM2_COLS = [
    "Completeness",
    "Contamination",
    "Completeness_Model_Used",
    "Translation_Table_Used",
    "Coding_Density",
    "Contig_N50",
    "Average_Gene_Length",
    "Genome_Size",
    "GC_Content",
    "Total_Coding_Sequences",
    "Total_Contigs",
    "Max_Contig_Length",
]

NUMERIC_COLS = [c for c in CHECKM2_COLS if c not in
                ("Completeness_Model_Used", "Translation_Table_Used")]

# Columns kept in the output for downstream compatibility, left empty
PLACEHOLDER_COLS_PRE = [
    "Organism Name",
    "Organism Infraspecific Names Strain",
    "Organism Infraspecific Names Isolate",
    "Annotation Name",
    "Assembly Level",
    "Assembly Release Date",
    "WGS project accession",
    "Assembly Stats Number of Scaffolds",
]
PLACEHOLDER_COLS_POST = ["batch_date", "Isolate", "WGS Project Accession"]

GENOME_EXT_REGEX = r"\.(?:fa|fna|fasta|fas)(?:\.gz)?$"


# ---------------------------------------------------------
# Helpers
# ---------------------------------------------------------

def log_msg(*msg):
    with open(log_path, "a") as f:
        f.write(" ".join(map(str, msg)) + "\n")


def read_table(path, label, sep):
    """Read a delimited table as text (separator is explicit, not sniffed)."""
    df = pd.read_csv(
        path,
        sep=sep,
        dtype=str,             # read everything as text; coerce numerics later
        skipinitialspace=True,
        encoding_errors="replace",
    )
    if df.empty:
        raise ValueError(f"{label} table is empty: {path}")
    return df


def clean_table(df):
    """Strip whitespace from headers and cell values; blank strings -> NA."""
    df = df.copy()
    df.columns = df.columns.astype(str).str.strip()
    for col in df.columns:
        df[col] = df[col].astype("string").str.strip()
    df = df.replace({"": pd.NA})
    df = df.dropna(how="all")
    return df


def genome_key(series):
    """Normalised key for matching genome names (strip whitespace + fasta extension)."""
    return (
        series.astype("string")
        .str.strip()
        .str.replace(GENOME_EXT_REGEX, "", regex=True)
    )


def require_columns(df, cols, label):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(
            f"{label} is missing required column(s): {missing}. "
            f"Found: {list(df.columns)}"
        )


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():
    log_msg("Starting CheckM2 + taxonomy processing")
    log_msg(f"CheckM2 input: {checkm2_file}")
    log_msg(f"Taxonomy input: {taxonomy_file}")
    log_msg(f"Unfiltered output: {output_unfiltered}")
    log_msg(f"Filtered output: {output_filtered}")

    # ---- Read CheckM2 ----
    checkm2 = clean_table(read_table(checkm2_file, "CheckM2", sep="\t"))
    log_msg(f"Read CheckM2 table: {len(checkm2)} rows, {len(checkm2.columns)} columns")

    require_columns(checkm2, ["Assembly_Name"] + CHECKM2_COLS, "CheckM2 table")

    checkm2 = checkm2.dropna(subset=["Assembly_Name"])
    dup_checkm2 = checkm2["Assembly_Name"].duplicated().sum()
    if dup_checkm2:
        log_msg(f"WARNING: {dup_checkm2} duplicated Assembly_Name values in CheckM2 table; keeping first")
        checkm2 = checkm2.drop_duplicates(subset="Assembly_Name", keep="first")

    # ---- Read taxonomy ----
    taxonomy = clean_table(read_table(taxonomy_file, "Taxonomy", sep=","))
    log_msg(f"Read taxonomy table: {len(taxonomy)} rows, {len(taxonomy.columns)} columns")

    require_columns(taxonomy, ["genome_file"], "Taxonomy table")

    absent_ranks = [r for r in TAX_RANKS if r not in taxonomy.columns]
    if absent_ranks:
        log_msg(f"WARNING: taxonomy table lacks rank column(s) {absent_ranks}; filling with NA")
        for r in absent_ranks:
            taxonomy[r] = pd.NA

    # ---- Build matching keys ----
    checkm2["_key"] = genome_key(checkm2["Assembly_Name"])
    taxonomy["_key"] = genome_key(taxonomy["genome_file"])

    taxonomy = taxonomy.dropna(subset=["_key"])
    dup_tax = taxonomy["_key"].duplicated().sum()
    if dup_tax:
        log_msg(f"WARNING: {dup_tax} duplicated genome entries in taxonomy table; keeping first")
        taxonomy = taxonomy.drop_duplicates(subset="_key", keep="first")

    # ---- Merge ----
    log_msg("Merging CheckM2 results with taxonomy")
    n_before = len(checkm2)
    df = checkm2.merge(
        taxonomy[["_key"] + TAX_RANKS],
        on="_key",
        how="left",
        validate="m:1",
    )
    assert len(df) == n_before, "Merge changed the number of rows"
    log_msg(f"Merge complete: {len(df)} rows")

    missing_taxonomy = df["domain"].isna().sum()
    log_msg(f"Genomes with missing taxonomy: {missing_taxonomy}/{len(df)}")

    unused_tax = len(set(taxonomy["_key"]) - set(df["_key"]))
    if unused_tax:
        log_msg(f"NOTE: {unused_tax} taxonomy entries had no matching CheckM2 genome")

    # ---- Numeric coercion ----
    for col in NUMERIC_COLS:
        coerced = pd.to_numeric(df[col], errors="coerce")
        n_bad = coerced.isna().sum() - df[col].isna().sum()
        if n_bad:
            log_msg(f"WARNING: {n_bad} non-numeric values in '{col}' set to NA")
        df[col] = coerced

    # ---- Construct output table ----
    log_msg("Constructing output table")

    out = pd.DataFrame(index=df.index)
    for col in PLACEHOLDER_COLS_PRE:
        out[col] = ""
    out["genome_file"] = df["_key"]          # Assembly_Name (extension stripped), kept even without a taxonomy match
    for col in CHECKM2_COLS:
        out[col] = df[col]
    for col in PLACEHOLDER_COLS_POST:
        out[col] = ""
    for col in TAX_RANKS:
        out[col] = df[col]

    # ---- Write unfiltered ----
    log_msg(f"Writing unfiltered table: {output_unfiltered}")
    out.to_csv(output_unfiltered, index=False)
    log_msg(f"Wrote unfiltered table with {len(out)} genomes")

    # ---- Filter ----
    log_msg(
        f"Applying CheckM2 thresholds: completeness >= {min_completeness}, "
        f"contamination <= {max_contamination}"
    )

    n_nan_qc = (out["Completeness"].isna() | out["Contamination"].isna()).sum()
    if n_nan_qc:
        log_msg(f"WARNING: {n_nan_qc} genomes have missing completeness/contamination and will be excluded")

    filtered = out[
        (out["Completeness"] >= min_completeness)
        & (out["Contamination"] <= max_contamination)
    ].copy()

    log_msg(f"Writing filtered table: {output_filtered}")
    filtered.to_csv(output_filtered, index=False)
    log_msg(f"Filtered genomes: {len(filtered)}/{len(out)}")

    if filtered.empty:
        log_msg("WARNING: no genomes passed the filter")

    log_msg("CheckM2 + taxonomy processing completed successfully")


try:
    main()
except Exception:
    log_msg("ERROR: processing failed")
    log_msg(traceback.format_exc())
    raise