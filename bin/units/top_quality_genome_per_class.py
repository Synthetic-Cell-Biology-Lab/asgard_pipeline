import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome_file",        required=True)
    parser.add_argument("--taxa_levels",        required=True)
    parser.add_argument("--completeness_col",   default="Completeness")
    parser.add_argument("--contamination_col",  default="Contamination")
    parser.add_argument("--output",             required=True)
    args = parser.parse_args()

    taxa_levels      = [t.strip() for t in args.taxa_levels.split(",")]
    completeness_col = args.completeness_col
    contamination_col = args.contamination_col

    genome_meta = pd.read_csv(args.genome_file)

    # Validate columns exist
    for col in [completeness_col, contamination_col]:
        if col not in genome_meta.columns:
            raise ValueError(f"Column '{col}' not found in genome file. "
                             f"Available: {genome_meta.columns.tolist()}")

    records = []

    for taxa_level in taxa_levels:
        if taxa_level not in genome_meta.columns:
            print(f"Warning: taxa level '{taxa_level}' not found, skipping.")
            continue

        for cls in genome_meta[taxa_level].dropna().unique():
            subset = genome_meta[genome_meta[taxa_level] == cls].copy()

            # Sort by completeness desc, contamination asc — pick top
            best = (
                subset
                .sort_values(
                    [completeness_col, contamination_col],
                    ascending=[False, True]
                )
                .iloc[0]
            )

            records.append({
                "taxa_level":       taxa_level,
                "class":            cls,
                "genome_file":      best["genome_file"],
                completeness_col:   best[completeness_col],
                contamination_col:  best[contamination_col],
            })

    result = pd.DataFrame(records)
    result.to_csv(args.output, index=False)
    print(f"Written to {args.output}")

if __name__ == "__main__":
    main()