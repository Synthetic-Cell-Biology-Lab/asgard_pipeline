# scripts/hyphy_parser.py

import json
import pandas as pd
import argparse
from pathlib import Path


def parse_fel(json_file):
    data = json.load(open(json_file))
    rows = data["MLE"]["content"]["0"]

    df = pd.DataFrame(rows, columns=[
        "alpha", "beta", "alpha_se", "beta_se",
        "LRT", "pvalue"
    ])
    df["site"] = df.index + 1
    df["selection"] = df["beta"] - df["alpha"]
    return df


def parse_meme(json_file):
    data = json.load(open(json_file))
    rows = data["MLE"]["content"]["0"]

    df = pd.DataFrame(rows, columns=[
        "alpha", "beta_minus", "p_minus",
        "beta_plus", "p_plus",
        "LRT", "pvalue"
    ])
    df["site"] = df.index + 1
    df["episodic"] = df["pvalue"] < 0.05
    return df


def parse_fubar(json_file):
    data = json.load(open(json_file))
    rows = data["MLE"]["content"]["0"]

    df = pd.DataFrame(rows, columns=[
        "alpha", "beta",
        "post_neg", "post_pos"
    ])
    df["site"] = df.index + 1
    df["selection"] = df["beta"] - df["alpha"]
    df["selected"] = df["post_pos"] > 0.9
    return df


def parse_absrel(json_file):
    data = json.load(open(json_file))
    branches = data["branch attributes"]["0"]

    records = []
    for branch, vals in branches.items():
        records.append({
            "branch": branch,
            "pvalue": vals.get("Corrected P-value"),
            "significant": vals.get("Corrected P-value", 1) < 0.05
        })

    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fel")
    parser.add_argument("--meme")
    parser.add_argument("--fubar")
    parser.add_argument("--absrel")
    parser.add_argument("--outdir")
    parser.add_argument("--prefix")

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.fel:
        df = parse_fel(args.fel)
        df.to_csv(outdir / f"{args.prefix}_FEL.csv", index=False)

    if args.meme:
        df = parse_meme(args.meme)
        df.to_csv(outdir / f"{args.prefix}_MEME.csv", index=False)

    if args.fubar:
        df = parse_fubar(args.fubar)
        df.to_csv(outdir / f"{args.prefix}_FUBAR.csv", index=False)

    if args.absrel:
        df = parse_absrel(args.absrel)
        df.to_csv(outdir / f"{args.prefix}_aBSREL.csv", index=False)


if __name__ == "__main__":
    main()