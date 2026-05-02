import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--synteny_csv",  required=True)
    parser.add_argument("--genome_file",  required=True)
    parser.add_argument("--taxa_levels",  required=True)
    parser.add_argument("--output",       required=True)
    args = parser.parse_args()

    taxa_levels = [t.strip() for t in args.taxa_levels.split(",")]

    synteny     = pd.read_csv(args.synteny_csv)
    genome_meta = pd.read_csv(args.genome_file)

    # Count copies of each Manual_annotation per genome
    counts = (
        synteny
        .groupby(["genome_file", "Manual_annotation"])
        .size()
        .reset_index(name="copy_count")
    )

    # Merge taxa metadata
    counts = counts.merge(genome_meta, on="genome_file", how="left")

    records = []
    chosen_genomes = set()  # grows as we assign genomes

    for taxa_level in taxa_levels:
        for cls in counts[taxa_level].dropna().unique():
            subset = counts[counts[taxa_level] == cls].copy()

            for annotation in subset["Manual_annotation"].unique():
                ann_subset = subset[subset["Manual_annotation"] == annotation].copy()

                max_copies = ann_subset["copy_count"].max()

                # Candidates tied at the max
                top_candidates = ann_subset[ann_subset["copy_count"] == max_copies]

                # Prefer already-chosen genome as tiebreaker
                already_chosen = top_candidates[
                    top_candidates["genome_file"].isin(chosen_genomes)
                ]

                if not already_chosen.empty:
                    winner = already_chosen.iloc[0]
                else:
                    winner = top_candidates.iloc[0]

                chosen_genomes.add(winner["genome_file"])

                records.append({
                    "taxa_level":        taxa_level,
                    "class":             cls,
                    "Manual_annotation": annotation,
                    "genome_file":       winner["genome_file"],
                    "copy_count":        winner["copy_count"],
                })

    result = pd.DataFrame(records)
    result = result[["taxa_level", "class", "Manual_annotation", "genome_file", "copy_count"]]
    result.to_csv(args.output, index=False)
    print(f"Written to {args.output}")

if __name__ == "__main__":
    main()