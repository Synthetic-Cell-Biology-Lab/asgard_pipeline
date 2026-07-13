from pathlib import Path
import pandas as pd


def collate_ffn(
    metadata_csv: str | Path,
    base_path: str | Path,
    output_fasta: str | Path,
):
    metadata_csv = Path(metadata_csv)
    base_path = Path(base_path)
    output_fasta = Path(output_fasta)

    df = pd.read_csv(metadata_csv)

    n_found = 0
    n_missing = 0

    with open(output_fasta, "w") as out_f:

        for genome in df["genome_file"]:

            genome = str(genome).strip()

            expected_ffn = base_path / genome / f"{genome}.ffn"

            if expected_ffn.exists():
                ffn_path = expected_ffn
            else:
                matches = list(base_path.rglob(f"{genome}.ffn"))

                if not matches:
                    print(f"Missing: {genome}")
                    n_missing += 1
                    continue

                ffn_path = matches[0]

            print(f"Adding: {ffn_path}")

            with open(ffn_path) as f:

                for line in f:

                    if line.startswith(">"):
                        header = line[1:].strip()
                        out_f.write(f">{header}\n")
                    else:
                        out_f.write(line)

            n_found += 1

    print(
        f"\nFinished.\n"
        f"Genomes added: {n_found}\n"
        f"Missing: {n_missing}\n"
        f"Output: {output_fasta}"
    )


if __name__ == "__main__":

    collate_ffn(
        metadata_csv="/home/anirudh/asgard_pipeline/database/collated/Version1/filtered/85comp10con/genome_file/jan2026_85comp10con_gf.csv",
        base_path="/home/anirudh/asgard_pipeline/database/cds_genomes",
        output_fasta="/home/anirudh/asgard_pipeline/database/collated/Version1/filtered/85comp10con/fna/v1_cp85_con10_dna.ffn",
    )
