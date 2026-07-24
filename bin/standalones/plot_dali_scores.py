"""
similarity_heatmap.py  –  publication-grade clustered heatmap

Usage:
    python similarity_heatmap.py --matrix MATRIX.tsv [--mapping MAPPING.tsv] [options]

Required arguments:
    --matrix PATH    Tab-separated file: first column is sample ID, remaining
                      columns are numeric similarity (or distance) values.
                      Column order must match row order (square matrix).
                      No header row.

Optional arguments:
    --mapping PATH         Tab-separated file: two columns with header
                            subject<TAB>accession. Maps each sample ID to a
                            group label. Samples sharing a label are averaged
                            before plotting. Omit (or pass --no-mapping) to
                            plot every row/column as-is.

    --out-prefix PREFIX    Stem for output files  [default: similarity_heatmap]
    --dpi DPI              Resolution for PNG     [default: 300]
    --figsize W H          Figure width height in inches  [default: 9 7.5]
    --title TEXT           Main title text
    --cbar-label TEXT      Colorbar label  [default: "Average similarity score"]
    --annotate             Print similarity values inside each cell
    --annot-decimals N     Decimal places for cell annotations  [default: 1]
    --mapping-order        Use the order accessions first appear in the mapping
                           file instead of clustering or matrix order
    --no-cluster           Skip hierarchical clustering, keep input order
    --no-mapping           Ignore mapping file / plot raw matrix
    --linkage METHOD       Linkage method for scipy  [default: average]
    --cmap NAME            Built-in matplotlib cmap OR "pub_teal" (default)
    --vmin FLOAT           Colormap minimum (auto if omitted)
    --vmax FLOAT           Colormap maximum (auto if omitted)
    --input-is-distance    Treat input matrix as a distance matrix (e.g. from
                            TM-align/US-align) rather than a similarity matrix;
                            skips the 1-normalised-similarity conversion before
                            clustering, and flips the diagonal-normalisation
                            assumption used for similarity matrices.
    --triangle {full,upper,lower}
                            Show only the upper or lower triangle of the
                            (symmetric) matrix, or the full matrix [default: full]
    --max-labels N          If more than N labels would be drawn, thin them to
                            at most N (evenly spaced) instead of overlapping
                            illegibly  [default: 80]
    --protein-list FILE     File containing one accession per line. Only these
                            proteins will be plotted.
    --export-matrix         Also write the final (averaged/reordered/filtered)
                            matrix out as PREFIX.tsv, for downstream use
                            (e.g. as a Snakemake rule output)
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
    expected_len = None
    with open(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            ids.append(parts[0])
            values = parts[1:]
            if expected_len is None:
                expected_len = len(values)
            elif len(values) != expected_len:
                sys.exit(
                    f"Matrix line {line_no} has {len(values)} value columns, "
                    f"expected {expected_len} (based on the first row). "
                    f"Offending line starts with: {parts[0]!r}"
                )
            try:
                rows.append([float(x) for x in values])
            except ValueError as e:
                sys.exit(f"Matrix line {line_no}: could not parse a number ({e}).")
    arr = np.array(rows)
    if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
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
            if subj in mapping and mapping[subj] != acc:
                print(
                    f"Warning: duplicate mapping entry for '{subj}' "
                    f"('{mapping[subj]}' -> '{acc}'); using the later value.",
                    file=sys.stderr,
                )
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

    Vectorized via a one-hot group-membership matrix G (n_samples x n_groups):
        group_sums   = G.T @ arr @ G
        group_counts = counts_row[:, None] * counts_col[None, :]
        group_avg    = group_sums / group_counts
    This avoids the O(n^2) pure-Python double loop of the previous
    implementation, which mattered once replicate counts got large.
    """
    missing = [s for s in ids if s not in mapping]
    if missing:
        sys.exit(
            f"Sample(s) not found in mapping file: {missing[:5]}"
            + (" ..." if len(missing) > 5 else "")
        )

    if desired_order is not None:
        acc_order = desired_order
    else:
        acc_order, seen = [], set()
        for s in ids:
            a = mapping[s]
            if a not in seen:
                acc_order.append(a)
                seen.add(a)

    acc_index = {a: i for i, a in enumerate(acc_order)}
    n_groups = len(acc_order)
    n_samples = len(ids)

    # One-hot membership matrix: G[i, g] = 1 if sample i belongs to group g
    G = np.zeros((n_samples, n_groups))
    for i, s in enumerate(ids):
        g = acc_index[mapping[s]]
        G[i, g] = 1.0

    counts = G.sum(axis=0)  # replicate count per group
    group_sums = G.T @ arr @ G
    group_counts = np.outer(counts, counts)
    with np.errstate(divide="ignore", invalid="ignore"):
        avg = group_sums / group_counts
    avg = np.nan_to_num(avg, nan=0.0)

    return acc_order, avg


# ── Clustering ────────────────────────────────────────────────────────────────
def cluster(labels, arr, method="average", input_is_distance=False):
    """Return (linkage_matrix, reordered_labels, reordered_array)."""
    if input_is_distance:
        dist = arr.copy()
    else:
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


# ── Label thinning ────────────────────────────────────────────────────────────
def thinned_labels(labels, max_labels):
    """
    If len(labels) > max_labels, keep at most max_labels evenly-spaced labels
    and blank the rest, so tick text stays legible on large matrices instead
    of overlapping into an unreadable smear.
    """
    n = len(labels)
    if n <= max_labels or max_labels <= 0:
        return list(labels)
    keep_idx = set(np.linspace(0, n - 1, max_labels, dtype=int).tolist())
    return [lbl if i in keep_idx else "" for i, lbl in enumerate(labels)]


# ── Triangle masking ──────────────────────────────────────────────────────────
def apply_triangle_mask(arr, mode):
    """Return a masked array (upper/lower triangle hidden as NaN) or the
    original array unchanged for mode == 'full'."""
    if mode == "full":
        return arr
    masked = arr.astype(float).copy()
    n = masked.shape[0]
    if mode == "upper":
        # keep upper triangle including diagonal, blank below it
        mask = np.tril(np.ones((n, n), dtype=bool), k=-1)
    elif mode == "lower":
        # keep lower triangle including diagonal, blank above it
        mask = np.triu(np.ones((n, n), dtype=bool), k=1)
    else:
        sys.exit(f"Unknown --triangle mode: {mode}")
    masked[mask] = np.nan
    return masked


# ── Plot ──────────────────────────────────────────────────────────────────────
def plot(labels, arr, Z, args):
    n = len(labels)
    cmap = make_pub_teal() if args.cmap == "pub_teal" else plt.get_cmap(args.cmap)
    cmap = cmap.copy() if hasattr(cmap, "copy") else cmap
    try:
        cmap.set_bad(color="white")
    except AttributeError:
        pass

    finite_vals = arr[np.isfinite(arr)] if np.isnan(arr).any() else arr
    vmin = args.vmin if args.vmin is not None else finite_vals.min()
    vmax = args.vmax if args.vmax is not None else finite_vals.max()

    w, h = args.figsize
    fig = plt.figure(figsize=(w, h), dpi=args.dpi)
    fig.patch.set_facecolor("white")

    if Z is not None:
        gs = GridSpec(
            2,
            3,
            figure=fig,
            width_ratios=[0.8, 6.5, 0.35],
            height_ratios=[1.0, 7.5],
            hspace=0.01,
            wspace=0.01,
            left=0.05,
            right=0.85,
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
        ax_dl.invert_yaxis()
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

    display_labels = thinned_labels(labels, args.max_labels)

    font_lbl = dict(fontsize=max(3.5, min(6.5, 60 / n)), fontfamily="DejaVu Sans")
    ax_heat.set_xticks(range(n))
    ax_heat.set_yticks(range(n))
    ax_heat.set_xticklabels(
        display_labels, rotation=45, ha="right", va="top", **font_lbl
    )
    ax_heat.set_yticklabels(display_labels, **font_lbl)
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
        cmap_obj = (
            make_pub_teal() if args.cmap == "pub_teal" else plt.get_cmap(args.cmap)
        )
        for r in range(n):
            for c in range(n):
                val = arr[r, c]
                if not np.isfinite(val):
                    continue  # masked-out triangle cell
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
    # Evenly spaced ticks that work for both integer-scaled and fractional
    # (e.g. 0-1 identity) value ranges, unlike the previous ceil/floor scheme
    # which produced an empty tick set whenever vmax - vmin < 1.
    cb.set_ticks(np.linspace(vmin, vmax, 6))

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
        if Z is not None:
            basis = "distance" if args.input_is_distance else "normalised similarity"
            sub = f"Hierarchical clustering ({args.linkage} linkage, 1 − {basis})"
        else:
            sub = "Input order preserved"
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


# ── Matrix export ─────────────────────────────────────────────────────────────
def export_matrix_tsv(path, labels, arr):
    with open(path, "w") as fh:
        for label, row in zip(labels, arr):
            values = "\t".join("" if not np.isfinite(v) else f"{v:.6g}" for v in row)
            fh.write(f"{label}\t{values}\n")


# ── CLI ───────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--matrix",
        required=True,
        help="Square similarity (or distance) matrix TSV (no header)",
    )
    p.add_argument(
        "--mapping",
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
    p.add_argument(
        "--input-is-distance",
        action="store_true",
        help="Treat --matrix as a distance matrix (e.g. TM-align/US-align "
        "output) rather than a similarity matrix for clustering purposes",
    )
    p.add_argument(
        "--triangle",
        choices=["full", "upper", "lower"],
        default="full",
        help="Show only one triangle of the symmetric matrix (default: full)",
    )
    p.add_argument(
        "--max-labels",
        type=int,
        default=80,
        help="Thin tick labels to at most this many if the matrix is larger "
        "(default: 80; use 0 to disable thinning)",
    )
    p.add_argument(
        "--protein-list",
        help="File containing one accession per line. Only these proteins will be plotted.",
    )
    p.add_argument(
        "--export-matrix",
        action="store_true",
        help="Also write the final matrix as PREFIX.tsv for downstream use",
    )
    args = p.parse_args()

    ids, arr = parse_matrix(args.matrix)

    if args.no_mapping or args.mapping is None:
        labels = ids
    else:
        mapping, acc_file_order = parse_mapping(args.mapping)
        desired = acc_file_order if args.mapping_order else None
        labels, arr = average_by_group(ids, arr, mapping, desired_order=desired)

    if args.protein_list:
        with open(args.protein_list) as fh:
            keep = {line.strip().split()[0] for line in fh if line.strip()}

        idx = [i for i, lbl in enumerate(labels) if lbl in keep]

        missing = keep - set(labels)
        if missing:
            print(
                f"Warning: {len(missing)} proteins from the list were not present "
                "in the matrix."
            )

        labels = [labels[i] for i in idx]
        arr = arr[np.ix_(idx, idx)]

    if args.mapping_order or args.no_cluster:
        Z, labels_final, arr_final = None, labels, arr
    else:
        Z, labels_final, arr_final = cluster(
            labels, arr, method=args.linkage, input_is_distance=args.input_is_distance
        )

    if args.export_matrix:
        export_matrix_tsv(f"{args.out_prefix}.tsv", labels_final, arr_final)
        print(f"Saved {args.out_prefix}.tsv")

    plot_arr = apply_triangle_mask(arr_final, args.triangle)
    fig = plot(labels_final, plot_arr, Z, args)

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
