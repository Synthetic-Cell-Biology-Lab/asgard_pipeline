"""
taxonomy_from_metadata  (Snakemake `script:` version)

GTT mode only. Select taxonomy ranks from GTDB metadata and restrict to
genomes surviving skder dereplication.

Expects the following from Snakemake:
    input.summary      genomes_summary.tsv
    input.derep_list   skder derep_genomes.txt
    output[0]          *_gtdbtk_classification_split.ar53.csv
    log[0]             log file
    wildcards.org
"""

import re
import time
import traceback
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------
# Snakemake objects
# ---------------------------------------------------------

summary_path = Path(snakemake.input.summary)
derep_path = Path(snakemake.input.derep_list)
out_path = Path(snakemake.output[0])
log_path = Path(snakemake.log[0])

org = snakemake.wildcards.org
threads = snakemake.threads


# ---------------------------------------------------------
# Constants
# ---------------------------------------------------------

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]

ACCESSION_CANDIDATES = ("Assembly Accession", "accession", "genome", "user_genome")
ASSEMBLY_NAME_COL = "ncbi_assembly_name"

FASTA_EXT_REGEX = r"\.(?:fa|fna|fasta|fas)(?:\.gz)?$"


# ---------------------------------------------------------
# Logging
# ---------------------------------------------------------

log_path.parent.mkdir(parents=True, exist_ok=True)
out_path.parent.mkdir(parents=True, exist_ok=True)
log_path.write_text("")  # start fresh


def logmsg(message):
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    with open(log_path, "a") as f:
        f.write(f"[{timestamp}] [taxonomy_from_metadata] {message}\n")


# ---------------------------------------------------------
# Helpers
# NOTE: pick_column / normalize_accession were defined in the Snakefile
# in the original rule. If your Snakefile versions do more than this
# (e.g. map GCF -> GCA), copy that logic here.
# ---------------------------------------------------------

def pick_column(df, candidates, label):
    for c in candidates:                       # exact match first
        if c in df.columns:
            return c

    lowered = {col.lower(): col for col in df.columns}
    for c in candidates:                       # then case-insensitive
        if c.lower() in lowered:
            return lowered[c.lower()]

    raise KeyError(
        f"Could not find a {label} column. Tried {list(candidates)}. "
        f"Columns found: {list(df.columns)}"
    )


def normalize_accession(value):
    """'GB_GCA_000000000.1' / 'RS_GCF_000000000.1' -> 'GCA_000000000.1' etc."""
    if pd.isna(value):
        return pd.NA
    v = re.sub(r"^(?:GB_|RS_)", "", str(value).strip())
    return v if v else pd.NA


def normalize_derep_id(line):
    """Accept bare IDs or full paths, with or without a FASTA extension."""
    name = Path(line.strip()).name
    return re.sub(FASTA_EXT_REGEX, "", name)


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():
    logmsg("=" * 80)
    logmsg("STARTING: taxonomy_from_metadata")
    logmsg(f"Wildcard org: {org}")
    logmsg(f"Working directory: {Path.cwd()}")
    logmsg(f"Threads allocated: {threads}")

    # ---- Inputs ----
    logmsg("-" * 80)
    logmsg("INPUT FILES")

    for label, p in [("summary", summary_path), ("derep_list", derep_path)]:
        logmsg(f"{label}: {p} (exists: {p.exists()})")
        if not p.exists():
            raise FileNotFoundError(f"{label} does not exist: {p}")
        stat = p.stat()
        logmsg(f"  Size: {stat.st_size:,} bytes; modified: {time.ctime(stat.st_mtime)}")

    # ---- Read metadata ----
    logmsg("-" * 80)
    logmsg("READING GTDB METADATA")

    df = pd.read_csv(summary_path, sep="\t", dtype=str)
    logmsg(f"Metadata shape: {df.shape}")

    # strip whitespace in headers and values; blank -> NA
    df.columns = df.columns.astype(str).str.strip()
    for col in df.columns:
        df[col] = df[col].str.strip()
    df = df.replace({"": pd.NA})

    if df.empty:
        raise ValueError(f"Metadata table is empty: {summary_path}")

    # ---- Accession column ----
    logmsg("-" * 80)
    logmsg("IDENTIFYING GENOME ACCESSION COLUMN")

    genome_col = pick_column(df, ACCESSION_CANDIDATES, "genome accession")
    logmsg(f"Selected accession column: {genome_col!r}")

    # ---- Validate ranks / required columns ----
    logmsg("-" * 80)
    logmsg("VALIDATING COLUMNS")

    missing = [c for c in RANKS + [ASSEMBLY_NAME_COL] if c not in df.columns]
    if missing:
        raise ValueError(
            f"GTDB metadata is missing required column(s): {missing}. "
            f"Columns found: {df.columns.tolist()}"
        )

    for rank in RANKS:
        logmsg(
            f"{rank}: unique={df[rank].nunique(dropna=True):,}, "
            f"missing={df[rank].isna().sum():,}"
        )

    # ---- Build user_genome ----
    logmsg("-" * 80)
    logmsg("BUILDING user_genome IDENTIFIERS")

    acc = df[genome_col].map(normalize_accession).astype("string")

    # NCBI file names replace whitespace in assembly names with underscores
    asm = (
        df[ASSEMBLY_NAME_COL]
        .astype("string")
        .str.replace(r"\s+", "_", regex=True)
    )

    n_no_asm = int(asm.isna().sum())
    if n_no_asm:
        logmsg(f"WARNING: {n_no_asm:,} rows have no '{ASSEMBLY_NAME_COL}'; "
               f"they cannot be matched and will be dropped")

    df["user_genome"] = acc + "_" + asm + "_genomic"   # NA if any part is NA

    logmsg("First 10 raw -> normalized:")
    for raw, norm in zip(df[genome_col].head(10), df["user_genome"].head(10)):
        logmsg(f"  {raw!r} -> {norm!r}")

    n_null = int(df["user_genome"].isna().sum())
    logmsg(f"Rows with null user_genome: {n_null:,}")
    df = df.dropna(subset=["user_genome"])

    n_dups = int(df["user_genome"].duplicated().sum())
    if n_dups:
        logmsg(f"WARNING: {n_dups:,} duplicated user_genome rows; keeping first")
        df = df.drop_duplicates(subset="user_genome", keep="first")

    logmsg(f"Unique usable identifiers: {df['user_genome'].nunique():,}")

    # ---- Derep list ----
    logmsg("-" * 80)
    logmsg("READING DEREPLICATION LIST")

    lines = [l for l in derep_path.read_text().splitlines() if l.strip()]
    derep = {normalize_derep_id(l) for l in lines}

    logmsg(f"Non-empty lines: {len(lines):,}; unique IDs: {len(derep):,}")
    for g in sorted(derep)[:10]:
        logmsg(f"  {g!r}")

    if not derep:
        raise ValueError(f"Dereplication list is empty: {derep_path}")

    # ---- Compare ----
    logmsg("-" * 80)
    logmsg("COMPARING METADATA AND DEREP IDENTIFIERS")

    metadata_ids = set(df["user_genome"])
    matched = metadata_ids & derep
    missing_from_metadata = derep - metadata_ids

    logmsg(f"Metadata IDs: {len(metadata_ids):,}; derep IDs: {len(derep):,}; "
           f"matched: {len(matched):,}")
    logmsg(f"Derep IDs absent from metadata: {len(missing_from_metadata):,}")
    logmsg(f"Metadata IDs not in derep list: {len(metadata_ids - derep):,}")

    for g in sorted(missing_from_metadata)[:20]:
        logmsg(f"  absent from metadata: {g!r}")

    # ---- Filter ----
    logmsg("-" * 80)
    logmsg("FILTERING METADATA TO DEREPLICATED GENOMES")

    rows_before = len(df)
    df = df[df["user_genome"].isin(derep)].copy()
    logmsg(f"Rows before: {rows_before:,}; after: {len(df):,}")

    if df.empty:
        raise ValueError(
            "No genomes remain after filtering: derep IDs do not match metadata IDs. "
            f"Example derep ID: {sorted(derep)[0]!r}; "
            f"example metadata ID: {sorted(metadata_ids)[0] if metadata_ids else None!r}"
        )

    # ---- Write ----
    logmsg("-" * 80)
    logmsg("WRITING OUTPUT")
    logmsg(f"Output path: {out_path}")

    df['genome_file'] = df['user_genome']

    df.to_csv(out_path, index=False)

    check = pd.read_csv(out_path, dtype=str)
    logmsg(f"Output size: {out_path.stat().st_size:,} bytes")
    logmsg(f"Output read-back shape: {check.shape}")

    if len(check) != len(df) or "user_genome" not in check.columns:
        raise ValueError("Output verification failed (row count or user_genome column).")

    logmsg("COMPLETED SUCCESSFULLY")
    logmsg("=" * 80)


try:
    main()
except Exception:
    logmsg("FAILED")
    logmsg(traceback.format_exc())
    raise