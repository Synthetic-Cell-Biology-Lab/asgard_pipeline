#!/usr/bin/env python3

import json
import pandas as pd
import argparse
from pathlib import Path
import numpy as np


# ---------------------- HELPERS ----------------------

def load_json(json_file):
    with open(json_file) as f:
        return json.load(f)


def parse_with_headers(json_file):
    data = load_json(json_file)

    try:
        content = data["MLE"]["content"]
        headers = [h[0] for h in data["MLE"]["headers"]]

        # --- Handle different formats ---
        if isinstance(content, dict):
            # FEL/MEME/FUBAR style
            rows = content.get("0", list(content.values())[0])
        elif isinstance(content, list):
            # SLAC style
            rows = content
        else:
            raise ValueError("Unknown content structure")

        diff = len(rows[0]) - len(headers)
        headers.extend([f"_{x}" for x in range(diff)])

    except Exception as e:
        raise ValueError(f"[ERROR] Unexpected format in {json_file}: {e}")

    df = pd.DataFrame(rows, columns=headers)
    df["site"] = df.index + 1

    return df


    
def find_column(df, keywords):
    for col in df.columns:
        for kw in keywords:
            if col.lower() == kw.lower():
                return col
    for col in df.columns:
        for kw in keywords:
            if kw.lower() in col.lower():
                return col
    return None


# ---------------------- PARSERS ----------------------

def parse_fel(json_file):
    df = parse_with_headers(json_file)

    alpha_col = find_column(df, ["alpha"])
    beta_col = find_column(df, ["beta"])
    p_col = find_column(df, ["p-value", "pval", "p"])

    if alpha_col and beta_col:
        df["omega"] = df[beta_col] / (df[alpha_col] + 1e-8)

        df["direction"] = np.where(
            df[beta_col] > df[alpha_col],
            "positive",
            "negative"
        )

    if p_col:
        pvals = df[p_col].replace(0, 1e-300)
        df["pvalue"] = pvals
        df["significant"] = pvals < 0.05
        df["neglog10_p"] = -np.log10(pvals)

    return df


def parse_slac(json_file):
    df = parse_with_headers(json_file)

    dn_col = find_column(df, ["dn", "dN"])
    ds_col = find_column(df, ["ds", "dS"])

    if dn_col and ds_col:
        df["omega"] = df[dn_col] / (df[ds_col] + 1e-8)

    # SLAC sometimes also has p-values
    p_col = find_column(df, ["p-value", "p"])

    if p_col:
        pvals = df[p_col].replace(0, 1e-300)
        df["pvalue"] = pvals
        df["significant"] = pvals < 0.05

    return df


def parse_meme(json_file):
    df = parse_with_headers(json_file)

    p_col = find_column(df, ["p-value"])

    if p_col:
        pvals = df[p_col].replace(0, 1e-300)
        df["pvalue"] = pvals
        df["neglog10_p"] = -np.log10(pvals)
        df["episodic"] = pvals < 0.05
    else:
        print(f"[WARN] No p-value column found in MEME: {json_file}")

    return df


def parse_fubar(json_file):
    df = parse_with_headers(json_file)

    pos_col = "Prob[alpha<beta]"
    neg_col = "Prob[alpha>beta]"

    if pos_col not in df.columns or neg_col not in df.columns:
        raise ValueError(f"[ERROR] Expected FUBAR columns not found in {json_file}")

    df["posterior_positive"] = df[pos_col]
    df["posterior_negative"] = df[neg_col]

    df["positive_selected"] = df[pos_col] > 0.9
    df["negative_selected"] = df[neg_col] > 0.9

    # Optional omega estimate
    alpha_col = find_column(df, ["alpha"])
    beta_col = find_column(df, ["beta"])

    if alpha_col and beta_col:
        df["omega"] = df[beta_col] / (df[alpha_col] + 1e-8)

    return df


def parse_absrel(json_file):
    data = load_json(json_file)

    try:
        branches = data["branch attributes"]["0"]
    except KeyError:
        raise ValueError(f"[ERROR] Unexpected aBSREL format: {json_file}")

    records = []

    for branch, vals in branches.items():
        pval = vals.get("Corrected P-value")

        records.append({
            "branch": branch,
            "pvalue": pval,
            "significant": (pval is not None and pval < 0.05),
            "rate_classes": vals.get("Rate classes")
        })

    return pd.DataFrame(records)


# ---------------------- MAIN ----------------------

def main():
    parser = argparse.ArgumentParser(
        description="Parse HyPhy JSON outputs into publication-ready CSV"
    )

    parser.add_argument("--fel")
    parser.add_argument("--meme")
    parser.add_argument("--fubar")
    parser.add_argument("--slac")
    parser.add_argument("--absrel")

    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", required=True)

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---------------- FEL ----------------
    if args.fel:
        print(f"[INFO] Parsing FEL")
        df = parse_fel(args.fel)
        df.to_csv(outdir / f"{args.prefix}_FEL.csv", index=False)

    # ---------------- SLAC ----------------
    if args.slac:
        print(f"[INFO] Parsing SLAC")
        df = parse_slac(args.slac)
        df.to_csv(outdir / f"{args.prefix}_SLAC.csv", index=False)

    # ---------------- MEME ----------------
    if args.meme:
        print(f"[INFO] Parsing MEME")
        df = parse_meme(args.meme)
        df.to_csv(outdir / f"{args.prefix}_MEME.csv", index=False)

    # ---------------- FUBAR ----------------
    if args.fubar:
        print(f"[INFO] Parsing FUBAR")
        df = parse_fubar(args.fubar)
        df.to_csv(outdir / f"{args.prefix}_FUBAR.csv", index=False)

    # ---------------- aBSREL ----------------
    if args.absrel:
        print(f"[INFO] Parsing aBSREL")
        df = parse_absrel(args.absrel)
        df.to_csv(outdir / f"{args.prefix}_aBSREL.csv", index=False)

    print("[INFO] Parsing complete.")


if __name__ == "__main__":
    main()