#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns


# ------------------ STYLE ------------------

sns.set_theme(style="whitegrid", context="talk")

PALETTE = {
    "FEL": "#1f77b4",        # blue
    "SLAC": "#2ca02c",       # green
    "MEME": "#d62728",       # red
    "FUBAR+": "#9467bd",     # purple
    "FUBAR-": "#8c564b"      # brown
}


# ------------------ HELPERS ------------------

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


# ------------------ PLOTTING ------------------

def plot_selection(fel, slac, meme, fubar, outpath, title):

    fig, ax = plt.subplots(figsize=(16,6))

    # ---------------- FEL ----------------
    if "omega" in fel.columns:
        fel_omega = fel["omega"].clip(0, 5)

        sns.lineplot(
            x=fel["site"],
            y=fel_omega,
            ax=ax,
            label="FEL",
            color=PALETTE["FEL"],
            linewidth=2
        )

    # ---------------- SLAC ----------------
    if slac is not None and "omega" in slac.columns:
        slac_omega = slac["omega"].clip(0, 5)
        smooth = pd.Series(slac_omega).rolling(7, center=True).mean()

        sns.lineplot(
            x=slac["site"],
            y=smooth,
            ax=ax,
            label="SLAC (smoothed)",
            color=PALETTE["SLAC"],
            linestyle="--",
            linewidth=2,
            alpha=0.9
        )

    # ---------------- MEME ----------------
    meme_p = find_column(meme, ["p"])
    if meme_p:
        meme_sig = meme[meme[meme_p] < 0.05]

        ax.vlines(
            meme_sig["site"],
            ymin=1,
            ymax=3,
            color=PALETTE["MEME"],
            alpha=0.6,
            linewidth=1.5,
            label="MEME (episodic)"
        )

    # ---------------- FUBAR ----------------
    pos_col = find_column(fubar, ["prob[alpha<beta]"])
    neg_col = find_column(fubar, ["prob[alpha>beta]"])

    if pos_col:
        pos_sites = fubar[fubar[pos_col] > 0.9]

        sns.scatterplot(
            x=pos_sites["site"],
            y=[2.6]*len(pos_sites),
            ax=ax,
            color=PALETTE["FUBAR+"],
            label="FUBAR +",
            s=40,
            edgecolor="black",
            linewidth=0.3
        )

    if neg_col:
        neg_sites = fubar[fubar[neg_col] > 0.9]

        sns.scatterplot(
            x=neg_sites["site"],
            y=[0.4]*len(neg_sites),
            ax=ax,
            color=PALETTE["FUBAR-"],
            label="FUBAR -",
            s=40,
            edgecolor="black",
            linewidth=0.3
        )

    # ---------------- Baseline ----------------
    ax.axhline(1, linestyle="--", color="black", alpha=0.7, linewidth=1.5)

    # ---------------- Labels ----------------
    ax.set_xlabel("Codon site", fontsize=14)
    ax.set_ylabel("dN/dS (ω)", fontsize=14)
    ax.set_title(title, fontsize=16, weight="bold")

    ax.set_ylim(0, 4)

    # Cleaner legend
    ax.legend(frameon=False, fontsize=11)

    sns.despine()

    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


# ------------------ MEME MANHATTAN ------------------

def plot_manhattan(meme, outpath, protein):

    p_col = find_column(meme, ["p-value", "pval", "p"])
    if not p_col:
        print("[WARN] No p-value column for MEME")
        return

    pvals = meme[p_col].replace(0, 1e-300)

    plt.figure(figsize=(16,5))

    sns.scatterplot(
        x=meme["site"],
        y=-np.log10(pvals),
        color="#444444",
        s=20
    )

    plt.axhline(-np.log10(0.05), linestyle="--", color="red", linewidth=1.5)

    plt.xlabel("Codon site")
    plt.ylabel("-log10(p-value)")
    plt.title(f"MEME significance: {protein}", weight="bold")

    sns.despine()
    plt.tight_layout()

    plt.savefig(outpath, dpi=300)
    plt.close()


# ------------------ SUMMARY ------------------

def generate_summary(fel, meme, fubar, absrel, outpath):

    fel_p = find_column(fel, ["p"])
    meme_p = find_column(meme, ["p"])
    fubar_pos = find_column(fubar, ["prob[alpha<beta]"])
    fubar_neg = find_column(fubar, ["prob[alpha>beta]"])

    summary = {}

    # FEL (use omega now)
    if "omega" in fel.columns and fel_p:
        sig = fel[fel[fel_p] < 0.05]
        summary["FEL_positive_sites"] = (sig["omega"] > 1).sum()
        summary["FEL_negative_sites"] = (sig["omega"] < 1).sum()

    # MEME
    if meme_p:
        summary["MEME_sites"] = (meme[meme_p] < 0.05).sum()

    # FUBAR
    if fubar_pos:
        summary["FUBAR_positive_sites"] = (fubar[fubar_pos] > 0.9).sum()

    if fubar_neg:
        summary["FUBAR_negative_sites"] = (fubar[fubar_neg] > 0.9).sum()

    # aBSREL
    if "significant" in absrel.columns:
        summary["aBSREL_branches"] = absrel["significant"].sum()

    pd.DataFrame([summary]).to_csv(outpath, index=False)


# ------------------ ITOL ------------------

def generate_itol(absrel, outpath):

    with open(outpath, "w") as f:
        f.write("DATASET_COLORSTRIP\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\taBSREL\n")
        f.write("COLOR\t#ff0000\n\n")
        f.write("DATA\n")

        for _, row in absrel.iterrows():
            if row.get("significant", False):
                f.write(f"{row['branch']}\t#ff0000\n")


# ------------------ MAIN ------------------

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--fel")
    parser.add_argument("--slac")
    parser.add_argument("--meme")
    parser.add_argument("--fubar")
    parser.add_argument("--absrel")
    parser.add_argument("--outdir")
    parser.add_argument("--prefix")

    args = parser.parse_args()

    outdir = Path(args.outdir)
    plot_dir = outdir / "plots"
    summary_dir = outdir / "summary"
    itol_dir = outdir / "itol"

    plot_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)
    itol_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    fel = pd.read_csv(args.fel)
    slac = pd.read_csv(args.slac) if args.slac else None
    meme = pd.read_csv(args.meme)
    fubar = pd.read_csv(args.fubar)
    absrel = pd.read_csv(args.absrel)

    # Plots
    plot_selection(
        fel, slac, meme, fubar,
        plot_dir / f"{args.prefix}_selection_plot.png",
        args.prefix
    )

    plot_manhattan(
        meme,
        plot_dir / f"{args.prefix}_manhattan.png",
        args.prefix
    )

    # Summary
    generate_summary(
        fel, meme, fubar, absrel,
        summary_dir / f"{args.prefix}_summary.csv"
    )

    # iTOL
    generate_itol(
        absrel,
        itol_dir / f"{args.prefix}_absrel.txt"
    )


if __name__ == "__main__":
    main()