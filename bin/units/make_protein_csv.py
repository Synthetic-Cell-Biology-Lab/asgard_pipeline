import pandas as pd
from pathlib import Path
from Bio import SeqIO


fasta = snakemake.input.fastas
genomes = snakemake.input.genomes
log = snakemake.log
output = snakemake.output


dfs = []


def log_msg(*msg):
    with open(log[0], "a") as f:
        f.write(" ".join(map(str, msg)) + "\n")


log_msg(f"Starting make_protein_csv for org: {snakemake.wildcards.org}")
log_msg(f"Genomes: {genomes}")
log_msg(f"Fastas: {fasta}")


for tsv_file, faa_file in zip(genomes, fasta):

    genome = Path(tsv_file).stem.replace(".tsv", "")

    log_msg("Entered the loop and processing genome: " + genome)

    # ---------------------------------------------------------
    # Parse FASTA once
    # ---------------------------------------------------------

    protein_info = {
        rec.id.split()[0]: {
            "protein_length": len(rec.seq),
            "protein_header": rec.description,
        }
        for rec in SeqIO.parse(faa_file, "fasta")
    }

    log_msg("Parsed FASTA for genome: " + genome)

    # ---------------------------------------------------------
    # Read annotation table
    # ---------------------------------------------------------

    with open(tsv_file) as f:
        first = next(f)

    if "bakta" in first.lower():

        df = pd.read_csv(
            tsv_file,
            sep="\t",
            skiprows=6,
            header=None
        )

        df.columns = [
            "Sequence Id",
            "Type",
            "Start",
            "Stop",
            "Strand",
            "Locus Tag",
            "Gene",
            "Product",
            "DbXrefs",
        ]

    else:
        df = pd.read_csv(
            tsv_file,
            sep="\t"
        )

    log_msg("Read annotation table for genome: " + genome)

    # ---------------------------------------------------------
    # Bakta
    # ---------------------------------------------------------

    if "Locus Tag" in df.columns:

        df = df.rename(
            columns={
                "Locus Tag": "locus_tag",
                "Type": "ftype",
                "Gene": "gene",
                "Product": "product",
            }
        )

        df["length_bp"] = (
            df["Stop"] - df["Start"] + 1
        ).abs()

        df["EC_number"] = ""
        df["COG"] = ""

    # ---------------------------------------------------------
    # Prokka
    # ---------------------------------------------------------

    elif "locus_tag" in df.columns:

        if "length_bp" not in df.columns:
            df["length_bp"] = (
                df["end"] - df["start"] + 1
            ).abs()

        if "ftype" not in df.columns:

            if "type" in df.columns:
                df["ftype"] = df["type"]
            else:
                df["ftype"] = ""

        if "gene" not in df.columns:
            df["gene"] = ""

        if "product" not in df.columns:
            df["product"] = ""

        if "EC_number" not in df.columns:
            df["EC_number"] = ""

        if "COG" not in df.columns:
            df["COG"] = ""

    else:

        raise ValueError(
            f"Could not determine annotation format for {tsv_file}\n"
            f"Columns found: {list(df.columns)}"
        )

    # ---------------------------------------------------------
    # Add protein information from FASTA
    # ---------------------------------------------------------

    df["protein_length"] = df["locus_tag"].map(
        lambda x: protein_info.get(
            x, {}
        ).get("protein_length")
    )

    df["protein_header"] = df["locus_tag"].map(
        lambda x: protein_info.get(
            x, {}
        ).get("protein_header")
    )

    # ---------------------------------------------------------
    # Construct output
    # ---------------------------------------------------------

    out = pd.DataFrame({
        "locus_tag": df["locus_tag"],
        "genome_file": genome,
        "protein_length": df["protein_length"],
        "protein_header": df["protein_header"],
        "ftype": df["ftype"],
        "length_bp": df["length_bp"],
        "gene": df["gene"],
        "EC_number": df["EC_number"],
        "COG": df["COG"],
        "product": df["product"],
        "predicted_structure": "",
        "IPS": "",
        "Manual_annotation": "",
    })

    dfs.append(out)


# ---------------------------------------------------------
# Combine all genomes
# ---------------------------------------------------------

if dfs:
    pd.concat(
        dfs,
        ignore_index=True
    ).to_csv(
        output[0],
        index=False
    )
else:
    raise ValueError("No genomes were processed.")