#!/usr/bin/env python3
"""
parse_hyphy.py  —  Parse HyPhy JSON outputs into publication-ready CSVs.

Supported methods: FEL, SLAC, MEME, FUBAR, aBSREL

Usage:
    python parse_hyphy.py --fel gene.FEL.json --meme gene.MEME.json \
                          --outdir results/ --prefix GENE1
"""

import json
import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(format="[%(levelname)s] %(message)s", level=logging.INFO)
log = logging.getLogger(__name__)

PVAL_THRESHOLD   = 0.05
PVAL_FLOOR       = 1e-300   # replace exact zeros before log-transform
OMEGA_EPSILON    = 1e-8     # avoid division by zero in dN/dS
FUBAR_THRESHOLD  = 0.9      # posterior probability cutoff


# ---------------------------------------------------------------------------
# JSON loading
# ---------------------------------------------------------------------------

def load_json(path: str | Path) -> dict:
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# MLE table extraction
#
# HyPhy uses at least three different content layouts depending on version
# and method. All are normalised here into a list-of-rows (row-major).
#
# Layout A — list of column-arrays (row-major columnar, most common in
#             newer FEL/MEME/FUBAR):
#   content = [ [alpha_s1, alpha_s2, ...],   ← column 0 values
#               [beta_s1,  beta_s2,  ...],   ← column 1 values
#               ... ]
#   → transpose to get rows
#
# Layout B — dict of rows (integer-keyed row-major):
#   content = {"0": [col0, col1, ...], "1": [...], ...}
#   → dict.values() are the rows, already in row-major order
#
# Layout C — nested dict (older SLAC):
#   content = {"0": {"0": [col0, col1, ...], "1": [...], ...}}
#   → unwrap outer "0" key, then same as Layout B
# ---------------------------------------------------------------------------

def _is_column_array_layout(content: list, n_headers: int) -> bool:
    """
    Return True when a list content block is columnar (Layout A):
    outer list = one entry per column, inner list = one entry per site.

    The reliable signal is: len(content) == n_headers (the declared column
    count), which means each element is a column-array, not a row.
    Fall back to the inner > outer heuristic when n_headers is unavailable.
    """
    if not content or not isinstance(content[0], list):
        return False
    if not all(isinstance(v, list) for v in content):
        return False
    inner_len = len(content[0])
    outer_len = len(content)
    if n_headers > 0:
        # Primary test: outer dimension matches declared header count
        if outer_len == n_headers and all(len(v) == inner_len for v in content):
            return True
        return False
    # Fallback when n_headers unknown: require inner strictly larger than outer
    return (inner_len > outer_len
            and all(len(v) == inner_len for v in content))


def _extract_mle_rows(content, n_headers: int = 0) -> list[list]:
    """
    Normalise any HyPhy MLE content block into a list of rows.

    Known layouts (checked in order):

    A  list of column-arrays  [ [col0_s1,col0_s2,...], [col1_s1,...], ... ]
       → detected when len(content)==n_headers; transposed into row-major.

    B  list of rows            [ [col0,col1,...], [col0,col1,...], ... ]
       → returned as-is.

    C  dict of rows            {"0":[col0,col1,...], "1":[col0,col1,...], ...}
       → dict.values() returned as rows.

    D  dict of column-arrays   {"0":[col0_s1,...], "1":[col0_s2,...], ...}
       → values transposed into rows.

    E  SLAC/MEME partition     {"0": {"by-site": {"AVERAGED": [rows]}, ...}}
       → drills into partition → by-site → AVERAGED to extract rows.

    F  single-partition wrap   {"0": <Layout A or B list>}
       → unwraps outer dict and recurses on content["0"].
    """
    # ---- Layout A / B: plain list ------------------------------------------
    if isinstance(content, list):
        if _is_column_array_layout(content, n_headers):
            return [list(row) for row in zip(*content)]
        return content

    if not isinstance(content, dict):
        raise ValueError(f"Unrecognised MLE content type: {type(content)}")

    # ---- Layout E: partition dict with by-site/AVERAGED --------------------
    # Matches SLAC and newer MEME: {"0": {"by-site": {"AVERAGED": [...]}, ...}}
    for part_val in content.values():
        if isinstance(part_val, dict) and "by-site" in part_val:
            by_site = part_val["by-site"]
            if isinstance(by_site, dict) and "AVERAGED" in by_site:
                averaged = by_site["AVERAGED"]
                if isinstance(averaged, list) and averaged:
                    if _is_column_array_layout(averaged, n_headers):
                        return [list(row) for row in zip(*averaged)]
                    return averaged

    values = list(content.values())

    # ---- Layout F: single-partition wrap {"0": <list>} --------------------
    # Must be checked BEFORE C/D because {"0": [col_arrays]} has values[0]
    # as a list, which would otherwise be returned as a single row by Layout C.
    # Applies to MEME: {"0": [col_array0, col_array1, ...]}
    if len(content) == 1:
        inner = next(iter(values))
        if isinstance(inner, list):
            return _extract_mle_rows(inner, n_headers)

    # ---- Layout C / D: dict whose values are lists -------------------------
    if all(isinstance(v, list) for v in values):
        if _is_column_array_layout(values, n_headers):
            # Layout D: dict values are column-arrays → transpose
            return [list(row) for row in zip(*values)]
        # Layout C: dict values are rows
        return values

    raise ValueError(f"Unrecognised MLE content structure: {type(content)}")


def extract_mle_dataframe(data: dict, source: str) -> pd.DataFrame:
    """
    Pull MLE headers + content out of a HyPhy result dict and return a
    DataFrame with a 1-based 'site' column added.
    """
    try:
        mle     = data["MLE"]
        headers = [h[0] for h in mle["headers"]]
        rows    = _extract_mle_rows(mle["content"], n_headers=len(headers))
    except KeyError as exc:
        raise ValueError(f"Missing expected MLE key in {source}: {exc}") from exc

    if not rows:
        raise ValueError(f"No rows found in MLE content of {source}")

    n_declared = len(headers)
    n_cols     = len(rows[0])

    # If columns far exceed declared headers the layout was almost certainly
    # misdetected (e.g. columnar data returned un-transposed as rows).
    if n_cols > n_declared * 10:
        raise ValueError(
            f"{source}: row width ({n_cols}) is >10x declared header count "
            f"({n_declared}). MLE content layout was likely misdetected."
        )

    if n_cols > n_declared:
        extras = [f"extra_{i}" for i in range(n_cols - n_declared)]
        headers.extend(extras)
        log.warning("%s: %d undeclared column(s) padded (extra_0..extra_%d)",
                    source, len(extras), len(extras) - 1)
    elif n_cols < n_declared:
        headers = headers[:n_cols]

    df = pd.DataFrame(rows, columns=headers)
    df.insert(0, "site", range(1, len(df) + 1))
    return df


# ---------------------------------------------------------------------------
# Column search
# Exact match first, then substring — prevents accidental collisions.
# ---------------------------------------------------------------------------

def find_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lower_cols = {c.lower(): c for c in df.columns}
    for name in candidates:
        if name.lower() in lower_cols:
            return lower_cols[name.lower()]
    for name in candidates:
        for col in df.columns:
            if name.lower() in col.lower():
                return col
    return None


def require_column(df: pd.DataFrame, candidates: list[str], context: str) -> str:
    col = find_column(df, candidates)
    if col is None:
        raise ValueError(f"{context}: could not find column matching {candidates}. "
                         f"Available: {list(df.columns)}")
    return col


# ---------------------------------------------------------------------------
# Shared derived-column helpers
# ---------------------------------------------------------------------------

def add_omega(df: pd.DataFrame,
              alpha_candidates: list[str],
              beta_candidates:  list[str],
              context: str) -> pd.DataFrame:
    alpha_col = require_column(df, alpha_candidates, context)
    beta_col  = require_column(df, beta_candidates,  context)
    df["dS"]    = df[alpha_col]
    df["dN"]    = df[beta_col]
    df["omega"] = df[beta_col] / (df[alpha_col] + OMEGA_EPSILON)
    df["direction"] = np.where(df[beta_col] > df[alpha_col], "positive", "negative")
    return df


def add_pvalue_columns(df: pd.DataFrame,
                       candidates: list[str],
                       context: str,
                       episodic_label: str = "significant") -> pd.DataFrame:
    p_col = require_column(df, candidates, context)
    pvals = pd.to_numeric(df[p_col], errors="coerce").clip(lower=PVAL_FLOOR)
    df["pvalue"]       = pvals
    df["neglog10_p"]   = -np.log10(pvals)
    df[episodic_label] = pvals < PVAL_THRESHOLD
    return df


# ---------------------------------------------------------------------------
# Per-method parsers
# ---------------------------------------------------------------------------

def parse_fel(path: str | Path) -> pd.DataFrame:
    """
    Fixed Effects Likelihood — site-level dN/dS with p-values.
    Returns: site, dS, dN, omega, direction, pvalue, neglog10_p, significant
    """
    data = load_json(path)
    df   = extract_mle_dataframe(data, "FEL")
    df   = add_omega(df, ["alpha", "alpha;synonymous"],     ["beta", "beta;non-synonymous"], "FEL")
    df   = add_pvalue_columns(df, ["p-value", "pval", "p"], "FEL")
    return df[["site", "dS", "dN", "omega", "direction", "pvalue", "neglog10_p", "significant"]]


def parse_slac(path: str | Path) -> pd.DataFrame:
    """
    SLAC — site-level substitution counts and rates.
    Uses: MLE → headers + content → 0 → by-site → RESOLVED
    """
    data = load_json(path)

    try:
        mle = data["MLE"]
        headers = mle["headers"]
        content = mle["content"]["0"]["by-site"]["RESOLVED"]
    except KeyError as exc:
        raise ValueError(f"SLAC: missing expected key {exc} in {path}") from exc

    # Column names
    columns = [h[0] for h in headers]

    df = pd.DataFrame(content, columns=columns)

    # Normalize column names (important)
    df.columns = [c.strip() for c in df.columns]

    # Convert numerics
    for col in df.columns:
        try:
            df[col] = pd.to_numeric(df[col])
        except Exception:
            pass

    # --- CORRECT omega ---
    if "dS" in df.columns and "dN" in df.columns:
        df["omega"] = df["dN"] / df["dS"]
    else:
        raise ValueError("SLAC: dN/dS columns not found")

    # Rename for consistency
    df["dS"] = df.get("dS")
    df["dN"] = df.get("dN")

    # --- SLAC p-values ---
    pos_col = find_column(df, ["P [dN/dS > 1]"])
    neg_col = find_column(df, ["P [dN/dS < 1]"])

    if pos_col:
        df["pvalue_positive"] = pd.to_numeric(df[pos_col], errors="coerce")
    if neg_col:
        df["pvalue_negative"] = pd.to_numeric(df[neg_col], errors="coerce")

    # Optional significance (you choose how strict)
    if pos_col:
        df["positive_selected"] = df["pvalue_positive"] < PVAL_THRESHOLD
    if neg_col:
        df["negative_selected"] = df["pvalue_negative"] < PVAL_THRESHOLD


    df.insert(0, "site", range(1, len(df) + 1))

    return df[[
        "site" if "site" in df.columns else df.index.name or df.index,
        "dS", "dN", "omega",
        "pvalue_positive", "pvalue_negative",
        "positive_selected", "negative_selected"
    ]]


def parse_meme(path: str | Path) -> pd.DataFrame:
    """
    Mixed Effects Model of Evolution — episodic positive selection.
    Returns: site, pvalue, neglog10_p, episodic
    """
    data = load_json(path)
    df   = extract_mle_dataframe(data, "MEME")
    df   = add_pvalue_columns(df, ["p-value", "pval"], "MEME", episodic_label="episodic")
    return df[["site", "pvalue", "neglog10_p", "episodic"]]


def parse_fubar(path: str | Path) -> pd.DataFrame:
    """
    Fast Unconstrained Bayesian AppRoximation — posterior probabilities.
    Returns: site, dS, dN, omega, posterior_positive, posterior_negative,
             positive_selected, negative_selected
    """
    data = load_json(path)
    df   = extract_mle_dataframe(data, "FUBAR")

    pos_col = require_column(df, ["Prob[alpha<beta]", "prob[alpha<beta]"], "FUBAR")
    neg_col = require_column(df, ["Prob[alpha>beta]", "prob[alpha>beta]"], "FUBAR")

    df["posterior_positive"]  = pd.to_numeric(df[pos_col], errors="coerce")
    df["posterior_negative"]  = pd.to_numeric(df[neg_col], errors="coerce")
    df["positive_selected"]   = df["posterior_positive"] > FUBAR_THRESHOLD
    df["negative_selected"]   = df["posterior_negative"] > FUBAR_THRESHOLD

    df = add_omega(df, ["alpha"], ["beta"], "FUBAR")

    return df[["site", "dS", "dN", "omega",
               "posterior_positive", "posterior_negative",
               "positive_selected",  "negative_selected"]]


def parse_absrel(path: str | Path) -> pd.DataFrame:
    """
    Adaptive Branch-Site Random Effects Likelihood — branch-level selection.
    Returns full scalar summary per branch.
    """
    data = load_json(path)

    try:
        branches = data["branch attributes"]["0"]
    except KeyError as exc:
        raise ValueError(f"aBSREL: missing expected key {exc} in {path}") from exc

    records = []

    for branch, attrs in branches.items():

        corrected_p = attrs.get("Corrected P-value", 1.0)
        uncorrected_p = attrs.get("Uncorrected P-value", 1.0)

        record = {
            "branch": branch,

            # Core stats
            "pvalue": corrected_p,
            "uncorrected_pvalue": uncorrected_p,

            "significant": corrected_p < PVAL_THRESHOLD,
            "uncorrected_significant": uncorrected_p < PVAL_THRESHOLD,

            # Model info
            "rate_classes": attrs.get("Rate classes"),
            "LRT": attrs.get("LRT"),

            # Evolutionary parameters
            "baseline_omega": attrs.get("Baseline MG94xREV omega ratio"),
            "baseline_length": attrs.get("Baseline MG94xREV"),

            "full_model_length": attrs.get("Full adaptive model"),
            "nonsyn_subs": attrs.get("Full adaptive model (non-synonymous subs/site)"),
            "syn_subs": attrs.get("Full adaptive model (synonymous subs/site)"),
        }

        records.append(record)

    df = pd.DataFrame(records)

    return df

# ---------------------------------------------------------------------------
# Registry — maps CLI flag name to parser function
# Adding a new method only requires one entry here.
# ---------------------------------------------------------------------------

PARSERS: dict[str, callable] = {
    "fel":    parse_fel,
    "slac":   parse_slac,
    "meme":   parse_meme,
    "fubar":  parse_fubar,
    "absrel": parse_absrel,
}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    for method in PARSERS:
        p.add_argument(f"--{method}", metavar="JSON",
                       help=f"Path to HyPhy {method.upper()} JSON output")
    p.add_argument("--outdir",  required=True, help="Output directory")
    p.add_argument("--prefix",  required=True, help="Filename prefix for CSVs")
    return p


def main():
    args   = build_parser().parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ran_any = False

    for method, parse_fn in PARSERS.items():
        path = getattr(args, method)
        if path is None:
            continue

        ran_any = True
        log.info("Parsing %s: %s", method.upper(), path)

        try:
            df = parse_fn(path)
        except Exception as exc:
            log.error("%s parsing failed: %s", method.upper(), exc)
            continue

        out_path = outdir / f"{args.prefix}_{method.upper()}.csv"
        df.to_csv(out_path, index=False)
        log.info("  → %s  (%d sites/branches)", out_path.name, len(df))

    if not ran_any:
        log.warning("No input files provided. Use --fel, --slac, --meme, --fubar, or --absrel.")


if __name__ == "__main__":
    main()