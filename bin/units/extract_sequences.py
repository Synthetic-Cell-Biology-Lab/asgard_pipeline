import os
import pandas as pd
from Bio import SeqIO
from pathlib import Path


def find_fasta_files(genome_dir):

    genome_dir = Path(genome_dir)

    faa_files = list(genome_dir.glob("*.faa"))
    fnn_files = list(genome_dir.glob("*.ffn"))

    if not faa_files:
        raise FileNotFoundError(f"No .faa file found in {genome_dir}")

    if not fnn_files:
        raise FileNotFoundError(f"No .fnn file found in {genome_dir}")

    # take the first match
    return str(faa_files[0]), str(fnn_files[0])


def load_existing_index(fasta):

    idx = fasta + ".idx"

    if not os.path.exists(idx):
        raise FileNotFoundError(f"Index file missing: {idx}")

    # this will load the existing SQLite index
    return SeqIO.index_db(idx, fasta, "fasta")


df = pd.read_csv(snakemake.input.subset)

cds_path = snakemake.params.cds_path

with open(snakemake.output.faa, "w") as faa_out, open(snakemake.output.ffn, "w") as fnn_out:

    grouped = df.groupby("genome_file")["locus_tag"].apply(list)

    for genome, loci in grouped.items():

        genome_dir = os.path.join(cds_path, genome)

        faa, ffn = find_fasta_files(genome_dir)

        if faa is None or ffn is None:
            print(f"Skipping {genome_dir}: missing FASTA files")
            continue

        faa_index = load_existing_index(faa)
        ffn_index = load_existing_index(ffn)

        for locus in loci:

            if locus in faa_index:
                SeqIO.write(faa_index[locus], faa_out, "fasta")

            if locus in ffn_index:
                SeqIO.write(ffn_index[locus], fnn_out, "fasta")