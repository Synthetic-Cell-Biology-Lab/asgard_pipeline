import argparse
import os
import traceback

import pandas as pd


# ---------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------

parser = argparse.ArgumentParser(
    description="Update master GTDB-Tk classification table"
)

parser.add_argument(
    "--taxonomy",
    nargs="+",
    required=True,
    help="GTDB-Tk taxonomy CSV files"
)

parser.add_argument(
    "--output-master",
    required=True,
    help="Output master GTDB-Tk classification CSV"
)

parser.add_argument(
    "--log",
    required=True,
    help="Log file"
)

def str2bool(value):
    """Parse booleans so '--update-mode', '--update-mode True' and
    '--update-mode false' all work (Snakemake often passes a value)."""

    if isinstance(value, bool):
        return value

    v = str(value).strip().lower()

    if v in ("true", "t", "yes", "y", "1"):
        return True

    if v in ("false", "f", "no", "n", "0"):
        return False

    raise argparse.ArgumentTypeError(
        f"Boolean value expected for --update-mode, got '{value}'"
    )


parser.add_argument(
    "--update-mode",
    nargs="?",
    const=True,
    default=False,
    type=str2bool,
    help="Update existing master table instead of replacing it "
         "(flag alone = True; also accepts True/False)"
)

args = parser.parse_args()


taxonomy_files = args.taxonomy
output_master = args.output_master
log_file = args.log
update_mode = args.update_mode


# ---------------------------------------------------------
# Logger
# ---------------------------------------------------------

def log_msg(*msg):

    log_dir = os.path.dirname(log_file)

    if log_dir:
        os.makedirs(
            log_dir,
            exist_ok=True
        )

    with open(log_file, "a") as f:
        f.write(
            " ".join(map(str, msg)) + "\n"
        )


# ---------------------------------------------------------
# Genome identifier homogenization
# ---------------------------------------------------------

ID_COLS = ["user_genome", "genome_file"]
GENOME_EXT_REGEX = r"\.(?:fa|fna|fasta|fas)(?:\.gz)?$"


def clean_genome_id(series):
    """Strip whitespace and the .fna (and .fa/.fasta/.gz variants) extension."""

    s = series.astype("string").str.strip()
    s = s.str.replace(GENOME_EXT_REGEX, "", regex=True)

    return s.replace("", pd.NA)


def homogenize_genome_ids(df, label):
    """
    Make 'user_genome' and 'genome_file' identical, clean, and extension-free.

    - Only one of the two present -> the other is created as a copy.
    - Both present -> both are cleaned; missing values in one are filled
      from the other; if they disagree, 'user_genome' wins.
    """

    df = df.copy()
    df.columns = df.columns.astype(str).str.strip()

    present = [c for c in ID_COLS if c in df.columns]

    if not present:

        log_msg(
            f"ERROR: Neither 'user_genome' nor 'genome_file' found in {label}"
        )

        raise ValueError(
            f"Could not determine genome identifier column in {label}. "
            f"Expected 'user_genome' and/or 'genome_file'. "
            f"Columns found: {list(df.columns)}"
        )

    for c in present:
        df[c] = clean_genome_id(df[c])

    if len(present) == 1:

        src = present[0]
        dst = "genome_file" if src == "user_genome" else "user_genome"

        log_msg(
            f"  {label}: only '{src}' present; creating '{dst}' as a copy"
        )

        df.insert(df.columns.get_loc(src) + 1, dst, df[src])

    else:

        both = df["user_genome"].notna() & df["genome_file"].notna()

        conflicts = (
            both & (df["user_genome"] != df["genome_file"]).fillna(False)
        ).sum()

        if conflicts:
            log_msg(
                f"  WARNING: {label}: {conflicts} rows where 'user_genome' "
                f"and 'genome_file' differ after cleaning; using 'user_genome'"
            )

        merged = df["user_genome"].fillna(df["genome_file"])

        df["user_genome"] = merged
        df["genome_file"] = merged

    log_msg(
        f"  {label}: genome identifiers homogenized "
        f"('user_genome' == 'genome_file', extensions stripped)"
    )

    return df


# ---------------------------------------------------------
# Start
# ---------------------------------------------------------

log_msg("=" * 70)
log_msg("Starting update_master_gtdbtk")
log_msg(f"Output: {output_master}")
log_msg(f"UPDATE_MODE: {update_mode}")
log_msg(
    f"Number of taxonomy input files: "
    f"{len(taxonomy_files)}"
)


try:

    # -----------------------------------------------------
    # Check input files
    # -----------------------------------------------------

    log_msg("Input taxonomy files:")

    for f in taxonomy_files:

        log_msg(f"  {f}")

        if not os.path.exists(f):

            log_msg(
                f"ERROR: Input file does not exist: {f}"
            )

            raise FileNotFoundError(
                f"GTDB-Tk taxonomy file not found: {f}"
            )

        log_msg(
            f"    Size: {os.path.getsize(f)} bytes"
        )


    # -----------------------------------------------------
    # Read input files
    # -----------------------------------------------------

    log_msg("Reading GTDB-Tk taxonomy files")

    dfs = []

    for f in taxonomy_files:

        log_msg(f"Reading: {f}")

        df = pd.read_csv(f)

        log_msg(
            f"  Rows: {len(df)}"
        )

        log_msg(
            f"  Columns: {list(df.columns)}"
        )


        # -------------------------------------------------
        # Homogenize genome identifier columns
        # -------------------------------------------------

        df = homogenize_genome_ids(df, label=f)


        log_msg(
            f"  Unique genomes in file: "
            f"{df['user_genome'].nunique()}"
        )

        dfs.append(df)


    # -----------------------------------------------------
    # Combine new data
    # -----------------------------------------------------

    log_msg(
        "Concatenating new GTDB-Tk results"
    )

    if not dfs:

        raise ValueError(
            "No GTDB-Tk taxonomy dataframes were loaded."
        )

    new_df = pd.concat(
        dfs,
        ignore_index=True
    )

    log_msg(
        f"Combined new data: "
        f"{len(new_df)} rows"
    )

    log_msg(
        f"Unique user_genome values: "
        f"{new_df['user_genome'].nunique()}"
    )


    # -----------------------------------------------------
    # Check for duplicate genomes in new data
    # -----------------------------------------------------

    new_duplicates = new_df[
        "user_genome"
    ].duplicated(
        keep=False
    ).sum()

    log_msg(
        f"Duplicate rows within new data: "
        f"{new_duplicates}"
    )


    # -----------------------------------------------------
    # Create master directory
    # -----------------------------------------------------

    master_dir = os.path.dirname(
        output_master
    )

    if master_dir:

        log_msg(
            f"Creating/checking master directory: "
            f"{master_dir}"
        )

        os.makedirs(
            master_dir,
            exist_ok=True
        )

    log_msg(
        "Master directory ready"
    )


    # -----------------------------------------------------
    # Load existing master
    # -----------------------------------------------------

    if update_mode and os.path.exists(output_master):

        log_msg(
            "UPDATE_MODE enabled and existing master found: "
            f"{output_master}"
        )

        old_df = pd.read_csv(
            output_master
        )

        log_msg(
            f"Existing master contains "
            f"{len(old_df)} rows"
        )

        log_msg(
            f"Existing master columns: "
            f"{list(old_df.columns)}"
        )


        # -------------------------------------------------
        # Homogenize existing master identifiers
        # (so old '.fna' IDs match new extension-free IDs)
        # -------------------------------------------------

        old_df = homogenize_genome_ids(
            old_df,
            label="existing master"
        )


        log_msg(
            f"Existing master unique genomes: "
            f"{old_df['user_genome'].nunique()}"
        )


        # -------------------------------------------------
        # Combine old + new
        # -------------------------------------------------

        combined = pd.concat(
            [
                old_df,
                new_df
            ],
            ignore_index=True
        )

        log_msg(
            f"Combined old + new data: "
            f"{len(combined)} rows"
        )

    else:

        if update_mode:

            log_msg(
                "UPDATE_MODE enabled, but no existing master "
                "was found. Creating a new master."
            )

        else:

            log_msg(
                "UPDATE_MODE disabled. Creating master "
                "from current taxonomy data only."
            )

        combined = new_df


    # -----------------------------------------------------
    # Remove duplicate genomes
    # -----------------------------------------------------

    before_dedup = len(combined)

    log_msg(
        f"Removing duplicate user_genome entries "
        f"(before: {before_dedup})"
    )

    duplicate_count = combined[
        "user_genome"
    ].duplicated(
        keep="last"
    ).sum()

    log_msg(
        f"Duplicate user_genome rows to remove: "
        f"{duplicate_count}"
    )

    combined = combined.drop_duplicates(
        subset="user_genome",
        keep="last"
    )

    log_msg(
        f"After deduplication: "
        f"{len(combined)} rows"
    )


    # -----------------------------------------------------
    # Final validation
    # -----------------------------------------------------

    if combined["user_genome"].isna().any():

        missing_count = combined[
            "user_genome"
        ].isna().sum()

        raise ValueError(
            f"Found {missing_count} rows with missing "
            "'user_genome' values after merging."
        )


    if combined["user_genome"].duplicated().any():

        raise ValueError(
            "Duplicate user_genome values remain after "
            "deduplication."
        )


    if not (combined["user_genome"] == combined["genome_file"]).all():

        raise ValueError(
            "'user_genome' and 'genome_file' are not identical "
            "after homogenization."
        )


    if combined["user_genome"].str.contains(
        GENOME_EXT_REGEX,
        regex=True
    ).any():

        raise ValueError(
            "Genome identifiers still contain a FASTA extension "
            "after homogenization."
        )


    # -----------------------------------------------------
    # Write master
    # -----------------------------------------------------

    log_msg(
        f"Writing master GTDB-Tk table: "
        f"{output_master}"
    )

    combined.to_csv(
        output_master,
        index=False
    )


    # -----------------------------------------------------
    # Verify output
    # -----------------------------------------------------

    if not os.path.exists(output_master):

        raise FileNotFoundError(
            f"Master output was not created: "
            f"{output_master}"
        )


    log_msg(
        f"Successfully wrote master file "
        f"({os.path.getsize(output_master)} bytes)"
    )


    # -----------------------------------------------------
    # Verify written file
    # -----------------------------------------------------

    log_msg(
        "Reading output back for verification"
    )

    verification_df = pd.read_csv(
        output_master
    )

    log_msg(
        f"Verified output rows: "
        f"{len(verification_df)}"
    )

    log_msg(
        f"Verified output columns: "
        f"{list(verification_df.columns)}"
    )


    for col in ID_COLS:

        if col not in verification_df.columns:

            raise ValueError(
                f"Output verification failed: "
                f"'{col}' column is missing."
            )


    if len(verification_df) != len(combined):

        raise ValueError(
            "Output verification failed: "
            f"expected {len(combined)} rows but found "
            f"{len(verification_df)} rows."
        )


    log_msg(
        f"Final number of genomes: "
        f"{verification_df['user_genome'].nunique()}"
    )


    # -----------------------------------------------------
    # Success
    # -----------------------------------------------------

    log_msg(
        "update_master_gtdbtk completed successfully"
    )

    log_msg("=" * 70)


except Exception as e:

    # -----------------------------------------------------
    # Error logging
    # -----------------------------------------------------

    log_msg("=" * 70)
    log_msg(
        "ERROR in update_master_gtdbtk"
    )

    log_msg(
        f"Exception type: "
        f"{type(e).__name__}"
    )

    log_msg(
        f"Exception message: {e}"
    )

    log_msg(
        "Traceback:"
    )

    log_msg(
        traceback.format_exc()
    )

    log_msg("=" * 70)

    raise