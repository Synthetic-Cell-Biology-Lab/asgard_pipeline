"""
similarity_heatmap.py  –  publication-grade clustered heatmap

Usage:
    python similarity_heatmap.py MATRIX.tsv MAPPING.tsv [options]

Positional arguments:
    matrix   Tab-separated file: first column is sample ID, remaining columns
             are numeric similarity values. Column order must match row order
             (square symmetric matrix). No header row.
    mapping  Tab-separated file: two columns with header  subject<TAB>accession
             Maps each sample ID to a group label. Samples sharing a label are
             averaged before plotting. Omit (or pass --no-mapping) to plot
             every row/column as-is.

Options:
    --out-prefix PREFIX   Stem for output files  [default: similarity_heatmap]
    --dpi DPI             Resolution for PNG     [default: 300]
    --figsize W H         Figure width height in inches  [default: 9 7.5]
    --title TEXT          Main title text
    --cbar-label TEXT     Colorbar label  [default: "Average similarity score"]
    --annotate            Print similarity values inside each cell
    --annot-decimals N    Decimal places for cell annotations  [default: 1]
    --mapping-order       Use the order accessions first appear in the mapping
                          file instead of clustering or matrix order
    --no-cluster          Skip hierarchical clustering, keep input order
    --no-mapping          Ignore mapping file / plot raw matrix
    --linkage METHOD      Linkage method for scipy  [default: average]
    --cmap NAME           Built-in matplotlib cmap OR "pub_teal" (default)
    --vmin FLOAT          Colormap minimum (auto if omitted)
    --vmax FLOAT          Colormap maximum (auto if omitted)
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform


# ── Colormap ──────────────────────────────────────────────────────────────────
def make_pub_teal():
    colors = [
        (0.98, 0.98, 1.00),
        (0.80, 0.90, 0.95),
        (0.45, 0.73, 0.82),
        (0.13, 0.53, 0.68),
        (0.02, 0.33, 0.50),
        (0.00, 0.17, 0.32),
    ]
    return LinearSegmentedColormap.from_list("pub_teal", colors, N=256)


# ── I/O helpers ───────────────────────────────────────────────────────────────
def parse_matrix(path):
    """Return (list_of_row_ids, 2-D float array)."""
    ids, rows = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            ids.append(parts[0])
            rows.append([float(x) for x in parts[1:]])
    arr = np.array(rows)
    if arr.shape[0] != arr.shape[1]:
        sys.exit(f"Matrix is {arr.shape[0]}×{arr.shape[1]} – must be square.")
    if len(ids) != arr.shape[0]:
        sys.exit("Row-ID count does not match matrix dimensions.")
    return ids, arr


def parse_mapping(path):
    """Return (dict {subject: accession}, list of accessions in file order)."""
    mapping = {}
    acc_file_order, seen = [], set()
    with open(path) as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                sys.exit(f"Mapping line {i+1} has fewer than 2 columns: {line!r}")
            subj, acc = parts[0], parts[1]
            if i == 0 and subj.lower() in ("subject", "sample", "id"):
                continue  # skip header
            mapping[subj] = acc
            if acc not in seen:
                acc_file_order.append(acc)
                seen.add(acc)
    return mapping, acc_file_order


# ── Averaging ─────────────────────────────────────────────────────────────────
def average_by_group(ids, arr, mapping, desired_order=None):
    """
    Collapse replicate rows/columns that share the same accession by averaging.
    desired_order: if given, use that accession order; otherwise use first-
                   appearance order in the matrix row IDs.
    Returns (ordered_accession_list, averaged_square_array).
    """
    if desired_order is not None:
        acc_order = desired_order
    else:
        acc_order, seen = [], set()
        for s in ids:
            if s not in mapping:
                sys.exit(f"Sample '{s}' not found in mapping file.")
            a = mapping[s]
            if a not in seen:
                acc_order.append(a)
                seen.add(a)

    for s in ids:
        if s not in mapping:
            sys.exit(f"Sample '{s}' not found in mapping file.")

    sums = defaultdict(lambda: defaultdict(float))
    counts = defaultdict(lambda: defaultdict(int))
    for ri, rs in enumerate(ids):
        ra = mapping[rs]
        for ci, cs in enumerate(ids):
            ca = mapping[cs]
            sums[ra][ca] += arr[ri, ci]
            counts[ra][ca] += 1

    n = len(acc_order)
    avg = np.array(
        [
            [
                sums[acc_order[r]][acc_order[c]] / counts[acc_order[r]][acc_order[c]]
                for c in range(n)
            ]
            for r in range(n)
        ]
    )
    return acc_order, avg


# ── Clustering ────────────────────────────────────────────────────────────────
def cluster(labels, arr, method="average"):
    """Return (linkage_matrix, reordered_labels, reordered_array)."""
    diag = np.diag(arr)
    # Normalise to [0,1] before converting to distance
    with np.errstate(divide="ignore", invalid="ignore"):
        sim_norm = arr / np.outer(np.sqrt(diag), np.sqrt(diag))
    sim_norm = np.nan_to_num(sim_norm, nan=0.0)
    dist = 1.0 - sim_norm
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2  # enforce symmetry
    dist = np.clip(dist, 0, None)
    Z = linkage(squareform(dist), method=method)
    order = leaves_list(Z)
    return Z, [labels[i] for i in order], arr[np.ix_(order, order)]


# ── Plot ──────────────────────────────────────────────────────────────────────
def plot(labels, arr, Z, args):
    n = len(labels)
    cmap = make_pub_teal() if args.cmap == "pub_teal" else plt.get_cmap(args.cmap)
    vmin = args.vmin if args.vmin is not None else arr.min()
    vmax = args.vmax if args.vmax is not None else arr.max()

    w, h = args.figsize
    fig = plt.figure(figsize=(w, h), dpi=args.dpi)
    fig.patch.set_facecolor("white")

    if Z is not None:
        gs = GridSpec(
            2,
            3,
            figure=fig,
            width_ratios=[1.0, 7.5, 0.35],
            height_ratios=[1.0, 7.5],
            hspace=0.01,
            wspace=0.015,
            left=0.02,
            right=0.92,
            top=0.93,
            bottom=0.20,
        )
        ax_dt = fig.add_subplot(gs[0, 1])
        ax_dl = fig.add_subplot(gs[1, 0])
        ax_heat = fig.add_subplot(gs[1, 1])
        ax_cb = fig.add_subplot(gs[1, 2])

        _dend_kwargs = dict(
            color_threshold=0,
            above_threshold_color="#444444",
            link_color_func=lambda k: "#444444",
        )
        dendrogram(Z, ax=ax_dt, orientation="top", no_labels=True, **_dend_kwargs)
        dendrogram(Z, ax=ax_dl, orientation="left", no_labels=True, **_dend_kwargs)
        ax_dt.set_axis_off()
        ax_dl.set_axis_off()
    else:
        gs = GridSpec(
            1,
            2,
            figure=fig,
            width_ratios=[7.5, 0.35],
            wspace=0.02,
            left=0.22,
            right=0.92,
            top=0.90,
            bottom=0.22,
        )
        ax_heat = fig.add_subplot(gs[0, 0])
        ax_cb = fig.add_subplot(gs[0, 1])

    im = ax_heat.imshow(
        arr, cmap=cmap, aspect="equal", vmin=vmin, vmax=vmax, interpolation="nearest"
    )

    font_lbl = dict(fontsize=max(3.5, min(6.5, 60 / n)), fontfamily="DejaVu Sans")
    ax_heat.set_xticks(range(n))
    ax_heat.set_yticks(range(n))
    ax_heat.set_xticklabels(labels, rotation=45, ha="right", va="top", **font_lbl)
    ax_heat.set_yticklabels(labels, **font_lbl)
    ax_heat.tick_params(axis="both", which="both", length=0, pad=3)

    for x in np.arange(-0.5, n, 1):
        ax_heat.axvline(x, color="white", lw=0.3)
        ax_heat.axhline(x, color="white", lw=0.3)
    for i in range(n):
        ax_heat.add_patch(
            plt.Rectangle(
                (i - 0.5, i - 0.5),
                1,
                1,
                fill=False,
                edgecolor="#222222",
                lw=0.5,
                zorder=3,
            )
        )
    ax_heat.spines[:].set_visible(False)

    # ── Cell annotations ──────────────────────────────────────────────────────
    if args.annotate:
        # Choose font size that fits: shrink aggressively for large matrices
        annot_fs = max(3.0, min(6.0, 52 / n))
        fmt = f".{args.annot_decimals}f"
        # Midpoint of colormap range for contrast flip
        mid = (vmin + vmax) / 2.0
        cmap_obj = (
            make_pub_teal() if args.cmap == "pub_teal" else plt.get_cmap(args.cmap)
        )
        # Sample colormap at midpoint to decide light/dark text threshold
        mid_rgb = cmap_obj(0.5)[:3]
        mid_lum = 0.2126 * mid_rgb[0] + 0.7152 * mid_rgb[1] + 0.0722 * mid_rgb[2]
        for r in range(n):
            for c in range(n):
                val = arr[r, c]
                t = (val - vmin) / (vmax - vmin) if vmax > vmin else 0.5
                cell_rgb = cmap_obj(t)[:3]
                lum = 0.2126 * cell_rgb[0] + 0.7152 * cell_rgb[1] + 0.0722 * cell_rgb[2]
                txt_color = "#1a1a1a" if lum > 0.45 else "#f0f0f0"
                ax_heat.text(
                    c,
                    r,
                    format(val, fmt),
                    ha="center",
                    va="center",
                    fontsize=annot_fs,
                    color=txt_color,
                    fontfamily="DejaVu Sans",
                )

    cb = fig.colorbar(im, cax=ax_cb, orientation="vertical")
    cb.set_label(args.cbar_label, fontsize=6, labelpad=6, fontfamily="DejaVu Sans")
    cb.ax.tick_params(labelsize=5.5, length=2, width=0.5)
    cb.outline.set_linewidth(0.5)
    cb.set_ticks(
        np.arange(np.ceil(vmin), np.floor(vmax) + 1, max(1, round((vmax - vmin) / 6)))
    )

    if args.title:
        fig.text(
            0.50,
            0.965,
            args.title,
            ha="center",
            va="top",
            fontsize=8,
            fontweight="bold",
            fontfamily="DejaVu Sans",
            color="#1a1a1a",
        )
        sub = (
            "Hierarchical clustering (average linkage, 1 − normalised similarity)"
            if Z is not None
            else "Input order preserved"
        )
        fig.text(
            0.50,
            0.945,
            sub,
            ha="center",
            va="top",
            fontsize=6,
            fontfamily="DejaVu Sans",
            color="#555555",
        )

    return fig


# ── CLI ───────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("matrix", help="Square similarity matrix TSV (no header)")
    p.add_argument(
        "mapping",
        nargs="?",
        default=None,
        help="subject→accession mapping TSV (header required)",
    )
    p.add_argument("--out-prefix", default="similarity_heatmap")
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument(
        "--figsize", type=float, nargs=2, default=[9, 7.5], metavar=("W", "H")
    )
    p.add_argument(
        "--title", default="Pairwise sequence similarity — averaged by accession"
    )
    p.add_argument("--cbar-label", default="Average similarity score")
    p.add_argument(
        "--annotate",
        action="store_true",
        help="Print similarity values inside each cell",
    )
    p.add_argument(
        "--annot-decimals",
        type=int,
        default=1,
        help="Decimal places for cell annotations (default: 1)",
    )
    p.add_argument(
        "--mapping-order",
        action="store_true",
        help="Order accessions as they appear in the mapping file",
    )
    p.add_argument("--no-cluster", action="store_true")
    p.add_argument("--no-mapping", action="store_true")
    p.add_argument(
        "--linkage",
        default="average",
        choices=[
            "average",
            "single",
            "complete",
            "ward",
            "weighted",
            "centroid",
            "median",
        ],
    )
    p.add_argument("--cmap", default="pub_teal")
    p.add_argument("--vmin", type=float, default=None)
    p.add_argument("--vmax", type=float, default=None)
    args = p.parse_args()

    ids, arr = parse_matrix(args.matrix)

    if args.no_mapping or args.mapping is None:
        labels = ids
        acc_file_order = None
    else:
        mapping, acc_file_order = parse_mapping(args.mapping)
        desired = acc_file_order if args.mapping_order else None
        labels, arr = average_by_group(ids, arr, mapping, desired_order=desired)

    if args.mapping_order or args.no_cluster:
        Z, labels_final, arr_final = None, labels, arr
    else:
        Z, labels_final, arr_final = cluster(labels, arr, method=args.linkage)

    fig = plot(labels_final, arr_final, Z, args)

    for ext in ("png", "pdf", "svg"):
        out = f"{args.out_prefix}.{ext}"
        kw = dict(bbox_inches="tight", facecolor="white")
        if ext == "png":
            kw["dpi"] = args.dpi
        fig.savefig(out, **kw)
        print(f"Saved {out}")
    plt.close()


if __name__ == "__main__":
    main()
