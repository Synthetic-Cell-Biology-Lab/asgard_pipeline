#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import seaborn as sns


# ------------------ STYLE ------------------

sns.set_theme(style="whitegrid", context="talk")
plt.rcParams["font.family"] = "DejaVu Sans"

PALETTE = {
    "FEL":    "#2166ac",   # blue
    "SLAC":   "#4dac26",   # green
    "MEME":   "#d7191c",   # red
    "FUBAR+": "#9467bd",   # purple
    "FUBAR-": "#8c564b",   # brown
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


def smooth(series, window=11):
    return pd.Series(series).rolling(window, center=True, min_periods=1).mean()


# ------------------ MAIN PLOT ------------------

def plot_selection(fel, slac, meme, fubar, outpath, title):
    """
    4-panel stacked layout — standard for HyPhy selection analyses:

      Panel 1 (top, tall)  : FEL ω line + SLAC ω smoothed line
      Panel 2              : FUBAR posterior probabilities (pos & neg)
      Panel 3              : MEME –log10(p) line with significance threshold
      Panel 4 (bottom,thin): Site-level rug / significance indicators
    """

    fig = plt.figure(figsize=(18, 14))
    gs  = gridspec.GridSpec(
        4, 1,
        height_ratios=[3, 2, 2, 0.6],
        hspace=0.08
    )

    ax_omega  = fig.add_subplot(gs[0])
    ax_fubar  = fig.add_subplot(gs[1], sharex=ax_omega)
    ax_meme   = fig.add_subplot(gs[2], sharex=ax_omega)
    ax_rug    = fig.add_subplot(gs[3], sharex=ax_omega)

    sites = fel["site"].values

    # ── Panel 1 : FEL & SLAC ω ────────────────────────────────────────────

    # FEL raw ω (clipped for readability)
    fel_omega = fel["omega"].clip(0, 5).values
    ax_omega.plot(sites, fel_omega,
                  color=PALETTE["FEL"], linewidth=1.5, alpha=0.85,
                  label="FEL ω")

    # shade FEL positive selection (ω > 1, significant)
    if "significant" in fel.columns and "omega" in fel.columns:
        sig_mask = (fel["significant"].values) & (fel["omega"].values > 1)
        ax_omega.fill_between(sites, 1, fel_omega,
                              where=sig_mask,
                              interpolate=True,
                              color=PALETTE["FEL"], alpha=0.18,
                              label="_nolegend_")

    # SLAC smoothed ω
    if slac is not None and "omega" in slac.columns:
        slac_omega = slac["omega"].clip(0, 5).values
        slac_sm    = smooth(slac_omega, window=15)
        ax_omega.plot(sites, slac_sm,
                      color=PALETTE["SLAC"], linewidth=2,
                      linestyle="--", alpha=0.9,
                      label="SLAC ω (smoothed)")

    ax_omega.axhline(1, color="black", linestyle=":", linewidth=1.4, alpha=0.7)
    ax_omega.set_ylim(0, 5)
    ax_omega.set_ylabel("dN/dS (ω)", fontsize=13)
    ax_omega.legend(frameon=False, fontsize=11, loc="upper right")
    ax_omega.set_title(title, fontsize=16, fontweight="bold", pad=12)
    ax_omega.tick_params(labelbottom=False)

    # ── Panel 2 : FUBAR posteriors ────────────────────────────────────────

    pos_col = find_column(fubar, ["posterior_positive", "prob[alpha<beta]"])
    neg_col = find_column(fubar, ["posterior_negative", "prob[alpha>beta]"])

    if pos_col:
        ax_fubar.plot(fubar["site"], fubar[pos_col],
                      color=PALETTE["FUBAR+"], linewidth=1.5,
                      label="FUBAR P(pos sel)")
        ax_fubar.fill_between(fubar["site"], 0, fubar[pos_col],
                              color=PALETTE["FUBAR+"], alpha=0.20)

    if neg_col:
        ax_fubar.plot(fubar["site"], fubar[neg_col],
                      color=PALETTE["FUBAR-"], linewidth=1.5,
                      label="FUBAR P(neg sel)")
        ax_fubar.fill_between(fubar["site"], 0, fubar[neg_col],
                              color=PALETTE["FUBAR-"], alpha=0.12)

    ax_fubar.axhline(0.9, color="grey", linestyle="--",
                     linewidth=1.2, alpha=0.7, label="0.9 threshold")
    ax_fubar.set_ylim(0, 1.05)
    ax_fubar.set_ylabel("Posterior prob.", fontsize=13)
    ax_fubar.legend(frameon=False, fontsize=10, loc="upper right")
    ax_fubar.tick_params(labelbottom=False)

    # ── Panel 3 : MEME –log10(p) ─────────────────────────────────────────

    p_col = find_column(meme, ["pvalue", "p-value", "p"])
    if p_col:
        pvals   = meme[p_col].replace(0, 1e-300)
        neg_log = -np.log10(pvals)

        ax_meme.plot(meme["site"], neg_log,
                     color=PALETTE["MEME"], linewidth=1.5, alpha=0.85,
                     label="MEME –log₁₀(p)")
        ax_meme.fill_between(meme["site"], 0, neg_log,
                             color=PALETTE["MEME"], alpha=0.15)

    thresh = -np.log10(0.05)
    ax_meme.axhline(thresh, color="black", linestyle="--",
                    linewidth=1.2, alpha=0.7, label="p = 0.05")
    ax_meme.set_ylabel("–log₁₀(p)", fontsize=13)
    ax_meme.set_ylim(0, None)
    ax_meme.legend(frameon=False, fontsize=10, loc="upper right")
    ax_meme.tick_params(labelbottom=False)

    # ── Panel 4 : significance rug ────────────────────────────────────────

    ax_rug.set_ylim(0, 3)
    ax_rug.set_yticks([])
    ax_rug.set_ylabel("", fontsize=11)

    # FEL positive
    if "significant" in fel.columns and "omega" in fel.columns:
        pos_sites = fel.loc[(fel["significant"]) & (fel["omega"] > 1), "site"]
        neg_sites = fel.loc[(fel["significant"]) & (fel["omega"] < 1), "site"]
        ax_rug.vlines(pos_sites, 2.1, 2.9,
                      color=PALETTE["FEL"], linewidth=1.2, alpha=0.7)
        ax_rug.vlines(neg_sites, 0.1, 0.9,
                      color="#aec7e8", linewidth=1.0, alpha=0.6)

    # MEME episodic
    if "episodic" in meme.columns:
        ep_sites = meme.loc[meme["episodic"], "site"]
        ax_rug.vlines(ep_sites, 1.1, 1.9,
                      color=PALETTE["MEME"], linewidth=1.2, alpha=0.7)

    # Legend patches
    rug_patches = [
        mpatches.Patch(color=PALETTE["FEL"],  label="FEL positive (sig)"),
        mpatches.Patch(color="#aec7e8",        label="FEL negative (sig)"),
        mpatches.Patch(color=PALETTE["MEME"],  label="MEME episodic"),
    ]
    ax_rug.legend(handles=rug_patches, frameon=False,
                  fontsize=9, loc="upper right", ncol=3)
    ax_rug.set_xlabel("Codon site", fontsize=13)

    # ── Final formatting ──────────────────────────────────────────────────

    for ax in [ax_omega, ax_fubar, ax_meme, ax_rug]:
        ax.set_xlim(sites.min(), sites.max())
        sns.despine(ax=ax, bottom=(ax != ax_rug))

    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[OK] Saved: {outpath}")


# ------------------ SUMMARY ------------------

def generate_summary(fel, meme, fubar, absrel, outpath):

    fel_p    = find_column(fel,  ["pvalue", "p"])
    meme_p   = find_column(meme, ["pvalue", "p"])
    fubar_pos = find_column(fubar, ["posterior_positive", "prob[alpha<beta]"])
    fubar_neg = find_column(fubar, ["posterior_negative", "prob[alpha>beta]"])

    summary = {}

    if "omega" in fel.columns and fel_p:
        sig = fel[fel[fel_p] < 0.05]
        summary["FEL_positive_sites"] = (sig["omega"] > 1).sum()
        summary["FEL_negative_sites"] = (sig["omega"] < 1).sum()

    if meme_p:
        summary["MEME_sites"] = (meme[meme_p] < 0.05).sum()

    if fubar_pos:
        summary["FUBAR_positive_sites"] = (fubar[fubar_pos] > 0.9).sum()
    if fubar_neg:
        summary["FUBAR_negative_sites"] = (fubar[fubar_neg] > 0.9).sum()

    if "significant" in absrel.columns:
        summary["aBSREL_branches"] = absrel["significant"].sum()

    pd.DataFrame([summary]).to_csv(outpath, index=False)
    print(f"[OK] Summary: {outpath}")


# ------------------ ITOL ------------------
def generate_itol(absrel, outpath):
    with open(outpath, "w") as f:
        # Header
        f.write("DATASET_STYLE\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\taBSREL\n")
        f.write("COLOR\t#ff0000\n\n")

        # Data section
        f.write("DATA\n")

        for _, row in absrel.iterrows():
            if row.get("significant", False):
                branch_id = row["branch"]

                # ID, TYPE, WHAT, COLOR, WIDTH, STYLE
                f.write(f"{branch_id}\tbranch\tnode\t#ff0000\t3\tnormal\n")

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

    outdir      = Path(args.outdir)
    plot_dir    = outdir / "plots"
    summary_dir = outdir / "summary"
    itol_dir    = outdir / "itol"

    for d in [plot_dir, summary_dir, itol_dir]:
        d.mkdir(parents=True, exist_ok=True)

    fel    = pd.read_csv(args.fel)
    slac   = pd.read_csv(args.slac) if args.slac else None
    meme   = pd.read_csv(args.meme)
    fubar  = pd.read_csv(args.fubar)
    absrel = pd.read_csv(args.absrel)

    plot_selection(
        fel, slac, meme, fubar,
        plot_dir / f"{args.prefix}_selection_plot.png",
        args.prefix
    )

    generate_summary(
        fel, meme, fubar, absrel,
        summary_dir / f"{args.prefix}_summary.csv"
    )

    generate_itol(
        absrel,
        itol_dir / f"{args.prefix}_absrel.txt"
    )


if __name__ == "__main__":
    main()