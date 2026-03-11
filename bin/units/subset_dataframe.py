#!/usr/bin/env python3

import sys
import pandas as pd

def main():

    if len(sys.argv) < 5:
        print("Usage: subset_dataframe.py <input_csv> <values> <column> <output_csv>")
        sys.exit(1)

    input_csv = sys.argv[1]
    values = sys.argv[2]
    column = sys.argv[3]
    output_csv = sys.argv[4]

    # Convert comma-separated string into list
    values = [v.strip() for v in values.split(",")]

    # Load dataframe
    df = pd.read_csv(input_csv)

    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in dataframe")

    # Subset dataframe
    subset = df[df[column].isin(values)]

    # Save result
    subset.to_csv(output_csv, index=False)

    print(f"Subset complete.")
    print(f"Input rows: {len(df)}")
    print(f"Output rows: {len(subset)}")
    print(f"Filtered values: {values}")

if __name__ == "__main__":
    main()