#!/usr/bin/env python3

import argparse
import math
import re
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

from Bio import SeqIO

AA = "ACDEFGHIKLMNPQRSTVWY"

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "database" / "protein_sets"


# ============================================================
# Utility functions
# ============================================================


def shannon_entropy(seq):
    """Calculate Shannon entropy of amino-acid composition."""

    if not seq:
        return 0.0

    counts = Counter(seq)
    n = len(seq)

    entropy = 0.0

    for count in counts.values():
        p = count / n
        entropy -= p * math.log2(p)

    return entropy


def longest_selected_run(seq, residues):
    """
    Find the longest consecutive run consisting only of
    the selected residues.

    `residues` can be either a string or a set.

    Examples:
        residues = "P"
        residues = "PG"
        residues = {"P", "G"}
    """

    if isinstance(residues, str):
        residues = set(residues.upper())
    else:
        residues = set(residues)

    best_length = 0
    best_start = None
    best_end = None

    current_length = 0
    current_start = None

    for i, aa in enumerate(seq):

        if aa in residues:

            if current_length == 0:
                current_start = i

            current_length += 1

            if current_length > best_length:
                best_length = current_length
                best_start = current_start
                best_end = i

        else:
            current_length = 0
            current_start = None

    if best_length == 0:
        return {
            "length": 0,
            "start": "",
            "end": "",
            "sequence": "",
        }

    return {
        "length": best_length,
        "start": best_start + 1,
        "end": best_end + 1,
        "sequence": seq[best_start : best_end + 1],
    }


def mask_sequence(seq, residues):
    """
    Replace all residues not in `residues` with X.

    Example:

        sequence = APPEPAPQPP
        residues = P

        output = XPPXPXPP

    For residues = PG:

        sequence = APPEGPAPQPP
        output   = XPPXGPXPP
    """

    return "".join(aa if aa in residues else "X" for aa in seq)


# ============================================================
# Sliding-window scanning
# ============================================================


def scan_windows(
    seq,
    window_size,
    step,
    min_fraction,
    residues,
):
    """
    Scan a sequence using a sliding window.

    `residues` defines which amino acids are considered.
    """

    hits = []

    residues = set(residues.upper())

    for start in range(
        0,
        len(seq) - window_size + 1,
        step,
    ):

        end = start + window_size

        window = seq[start:end]

        selected_count = sum(aa in residues for aa in window)

        fraction = selected_count / window_size

        if fraction < min_fraction:
            continue

        entropy = shannon_entropy(window)

        selected_run = longest_selected_run(
            window,
            residues,
        )

        hits.append(
            {
                "start": start + 1,
                "end": end,
                "length": window_size,
                "selected_count": selected_count,
                "fraction": fraction,
                "longest_selected_run": selected_run["length"],
                "entropy": entropy,
                "sequence": window,
                "masked_sequence": mask_sequence(
                    window,
                    residues,
                ),
            }
        )

    return hits


# ============================================================
# Merge windows
# ============================================================


def merge_hits(hits, max_gap=0):
    """
    Merge overlapping or closely adjacent windows.
    """

    if not hits:
        return []

    hits = sorted(
        hits,
        key=lambda x: (
            x["start"],
            x["end"],
        ),
    )

    merged = []

    current_start = hits[0]["start"]
    current_end = hits[0]["end"]

    for hit in hits[1:]:

        if hit["start"] <= current_end + max_gap + 1:

            current_end = max(
                current_end,
                hit["end"],
            )

        else:

            merged.append(
                {
                    "start": current_start,
                    "end": current_end,
                }
            )

            current_start = hit["start"]
            current_end = hit["end"]

    merged.append(
        {
            "start": current_start,
            "end": current_end,
        }
    )

    return merged


# ============================================================
# Region analysis
# ============================================================


def analyze_region(
    seq,
    start,
    end,
    global_fraction,
    residues,
):
    """Analyze a merged region."""

    residues = set(residues.upper())

    region = seq[start - 1 : end]

    length = len(region)

    selected_count = sum(aa in residues for aa in region)

    fraction = selected_count / length

    if global_fraction > 0:

        enrichment = fraction / global_fraction

    else:

        enrichment = float("inf")

    selected_run = longest_selected_run(
        region,
        residues,
    )

    if selected_run["start"] != "":

        absolute_run_start = start + selected_run["start"] - 1

        absolute_run_end = start + selected_run["end"] - 1

    else:

        absolute_run_start = ""
        absolute_run_end = ""

    return {
        "start": start,
        "end": end,
        "length": length,
        "selected_count": selected_count,
        "fraction": fraction,
        "enrichment": enrichment,
        "longest_selected_run": selected_run["length"],
        "longest_selected_start": absolute_run_start,
        "longest_selected_end": absolute_run_end,
        "longest_selected_sequence": selected_run["sequence"],
        "entropy": shannon_entropy(region),
        "sequence": region,
        "masked_sequence": mask_sequence(
            region,
            residues,
        ),
    }


# ============================================================
# Pattern search
# ============================================================


def pattern_to_regex(pattern):
    """
    Convert a protein pattern into a regular expression.

    Amino-acid letters and regex syntax are kept as-is; a
    wildcard `X` outside of a character class matches any
    standard amino acid.

    Examples:
        PXXPXXP                      ->  P[A-Z][A-Z]P[A-Z][A-Z]P
        PPP[PA]P((P[LGP])|([LG]P))   ->  unchanged
    """

    in_class = False

    regex_chars = []

    for char in pattern:

        if char == "[":

            in_class = True

        elif char == "]":

            in_class = False

        if char == "X" and not in_class:

            regex_chars.append("[A-Z]")

        else:

            regex_chars.append(char)

    return "".join(regex_chars)


def scan_pattern(seq, pattern_regex):
    """
    Find all non-overlapping pattern matches in a sequence.

    Returns a list of match dictionaries with 1-based start/end
    coordinates.
    """

    matches = []

    for match in pattern_regex.finditer(seq):

        matches.append(
            {
                "start": match.start() + 1,
                "end": match.end(),
                "length": match.end() - match.start(),
                "sequence": match.group(),
            }
        )

    return matches


# ============================================================
# Protein processing
# ============================================================


def process_protein(
    record,
    window_sizes,
    step,
    min_fraction,
    max_gap,
    residues,
    pattern_regex=None,
):
    """
    Process one protein.

    Returns:
        windows_by_size
        regions_by_size
        protein_summary
        pattern_matches
        pattern_hit_seqs
    """

    raw_seq = str(record.seq).upper()

    # Keep only standard amino acids.
    seq = "".join(aa for aa in raw_seq if aa in AA)

    if not seq:
        return {}, {}, None, [], []

    protein_length = len(seq)

    residue_set = set(residues.upper())

    global_selected_count = sum(aa in residue_set for aa in seq)

    global_fraction = global_selected_count / protein_length

    longest_selected = longest_selected_run(
        seq,
        residue_set,
    )

    windows_by_size = {}
    regions_by_size = {}

    for window_size in window_sizes:

        if len(seq) < window_size:
            continue

        hits = scan_windows(
            seq,
            window_size,
            step,
            min_fraction,
            residues,
        )

        for hit in hits:

            hit["protein_id"] = record.id

            hit["window_size"] = window_size

            hit["global_fraction"] = global_fraction

            if global_fraction > 0:

                hit["enrichment"] = hit["fraction"] / global_fraction

            else:

                hit["enrichment"] = float("inf")

        windows_by_size[window_size] = hits

        merged = merge_hits(
            hits,
            max_gap=max_gap,
        )

        regions = []

        for region in merged:

            result = analyze_region(
                seq,
                region["start"],
                region["end"],
                global_fraction,
                residues,
            )

            result["protein_id"] = record.id

            result["window_size"] = window_size

            result["global_fraction"] = global_fraction

            regions.append(result)

        regions_by_size[window_size] = regions

    pattern_matches = []

    pattern_hit_seqs = []

    if pattern_regex is not None:

        pattern_matches = scan_pattern(
            seq,
            pattern_regex,
        )

        if pattern_matches:

            pattern_hit_seqs.append(
                {
                    "protein_id": record.id,
                    "sequence": seq,
                }
            )

        for match in pattern_matches:

            match["protein_id"] = record.id

    protein_summary = {
        "protein_id": record.id,
        "length": protein_length,
        "total_selected": global_selected_count,
        "global_fraction": global_fraction,
        "longest_selected_run": longest_selected["length"],
        "longest_selected_start": longest_selected["start"],
        "longest_selected_end": longest_selected["end"],
        "longest_selected_sequence": longest_selected["sequence"],
    }

    return (
        windows_by_size,
        regions_by_size,
        protein_summary,
        pattern_matches,
        pattern_hit_seqs,
    )


# ============================================================
# TSV writer
# ============================================================


def write_tsv(
    path,
    rows,
    columns,
):
    """Write rows to TSV."""

    with open(path, "w") as out:

        out.write("\t".join(columns) + "\n")

        for row in rows:

            values = []

            for column in columns:

                value = row.get(
                    column,
                    "",
                )

                if isinstance(
                    value,
                    float,
                ):

                    if math.isinf(value):

                        value = "inf"

                    else:

                        value = f"{value:.6f}"

                values.append(str(value))

            out.write("\t".join(values) + "\n")


# ============================================================
# Directory creation
# ============================================================


def create_output_structure(
    output_dir,
    window_sizes,
):
    """Create organized output directories."""

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    rankings_dir = output_dir / "rankings"

    windows_dir = output_dir / "windows"

    regions_dir = output_dir / "regions"

    rankings_dir.mkdir(exist_ok=True)

    windows_dir.mkdir(exist_ok=True)

    regions_dir.mkdir(exist_ok=True)

    for window_size in window_sizes:

        (windows_dir / f"{window_size}aa").mkdir(exist_ok=True)

        (regions_dir / f"{window_size}aa").mkdir(exist_ok=True)


# ============================================================
# Log
# ============================================================


def write_log(
    path,
    args,
    protein_count,
    window_counts,
    region_counts,
    elapsed_seconds,
):
    """Write complete run information."""

    with open(path, "w") as log:

        log.write("============================================================\n")

        log.write("RESIDUE-RICH REGION SCANNER\n")

        log.write("============================================================\n\n")

        log.write(f"Run time: " f"{datetime.now().isoformat()}\n")

        log.write(f"Runtime: " f"{elapsed_seconds:.2f} seconds\n")

        log.write("\n")

        log.write("--------------------\n")

        log.write("INPUT\n")

        log.write("--------------------\n")

        log.write(f"Input FASTA: " f"{args.input}\n")

        log.write("\n")

        log.write("--------------------\n")

        log.write("OUTPUT\n")

        log.write("--------------------\n")

        log.write(f"Output directory: " f"{args.output_dir}\n")

        log.write("\n")

        log.write("--------------------\n")

        log.write("SCAN PARAMETERS\n")

        log.write("--------------------\n")

        log.write(f"Selected residues: " f"{args.residues}\n")

        if getattr(
            args,
            "pattern",
            None,
        ):

            log.write(
                f"Pattern: "
                f"{args.pattern.upper()}\n"
            )

        log.write(
            "Window sizes: "
            + ",".join(
                map(
                    str,
                    args.windows,
                )
            )
            + "\n"
        )

        log.write(f"Step size: " f"{args.step}\n")

        log.write(f"Minimum fraction: " f"{args.min_fraction}\n")

        log.write(f"Maximum merge gap: " f"{args.max_gap}\n")

        log.write("\n")

        log.write("--------------------\n")

        log.write("RUN SUMMARY\n")

        log.write("--------------------\n")

        log.write(f"Proteins scanned: " f"{protein_count}\n")

        log.write("\n")

        for window_size in args.windows:

            log.write(f"{window_size} aa windows:\n")

            log.write(
                f"    Selected-residue windows: "
                f"{window_counts.get(window_size, 0)}\n"
            )

            log.write(f"    Regions: " f"{region_counts.get(window_size, 0)}\n")

        log.write("\n")

        log.write("--------------------\n")

        log.write("COMMAND\n")

        log.write("--------------------\n")

        log.write(" ".join(sys.argv) + "\n")


# ============================================================
# Main
# ============================================================


def main():

    parser = argparse.ArgumentParser(
        description=("Scan protein FASTA files for " "residue-rich regions.")
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Protein FASTA file",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help=(
            "Output directory. "
            f"Default: {DEFAULT_OUTPUT_DIR}"
        ),
    )

    parser.add_argument(
        "-w",
        "--windows",
        nargs="+",
        type=int,
        default=[
            10,
            15,
            20,
            25,
            30,
            40,
            50,
        ],
        help="Window sizes",
    )

    parser.add_argument(
        "-s",
        "--step",
        type=int,
        default=1,
        help="Sliding-window step size",
    )

    parser.add_argument(
        "--min-fraction",
        type=float,
        default=0.30,
        help=("Minimum fraction of selected " "residues in a window"),
    )

    parser.add_argument(
        "--max-gap",
        type=int,
        default=0,
        help=("Maximum gap between windows " "when merging"),
    )

    parser.add_argument(
        "--residues",
        default="P",
        help=("Residues to consider. " "Examples: P, G, PAG. " "Default: P"),
    )

    parser.add_argument(
        "--pattern",
        default=None,
        help=(
            "Protein pattern to search for. "
            "Regex syntax is supported and X is a wildcard "
            "for any amino acid. "
            "Examples: PXXPXXP, "
            "PPP[PA]P((P[LGP])|([LG]P))"
        ),
    )

    args = parser.parse_args()

    # --------------------------------------------------------
    # Validate residues
    # --------------------------------------------------------

    residues = set(args.residues.upper())

    invalid_residues = residues - set(AA)

    if invalid_residues:

        parser.error(
            "Invalid amino-acid residue(s): " + ", ".join(sorted(invalid_residues))
        )

    if not residues:

        parser.error("--residues cannot be empty")

    # --------------------------------------------------------
    # Validate pattern
    # --------------------------------------------------------

    pattern = None

    pattern_regex = None

    if args.pattern:

        pattern = args.pattern.upper()

        try:

            pattern_regex = re.compile(
                pattern_to_regex(pattern)
            )

        except re.error as error:

            parser.error(
                f"Invalid pattern '{args.pattern}': {error}"
            )

    # --------------------------------------------------------
    # Validate parameters
    # --------------------------------------------------------

    if args.step < 1:

        parser.error("--step must be >= 1")

    if args.min_fraction < 0 or args.min_fraction > 1:

        parser.error("--min-fraction must be between 0 and 1")

    if args.max_gap < 0:

        parser.error("--max-gap must be >= 0")

    if any(window <= 0 for window in args.windows):

        parser.error("Window sizes must be > 0")

    # --------------------------------------------------------
    # Setup
    # --------------------------------------------------------

    start_time = datetime.now()

    output_dir = Path(args.output_dir)

    create_output_structure(
        output_dir,
        args.windows,
    )

    print("==================================================")

    print("Residue-Rich Region Scanner")

    print("==================================================")

    print(f"Input:          {args.input}")

    print(f"Output:         {output_dir}")

    print(f"Residues:       {args.residues.upper()}")

    if pattern:

        print(f"Pattern:        {pattern}")

    print(f"Windows:        {args.windows}")

    print(f"Step:            {args.step}")

    print(f"Minimum fraction: {args.min_fraction}")

    print()

    # --------------------------------------------------------
    # Storage
    # --------------------------------------------------------

    all_windows = {window_size: [] for window_size in args.windows}

    all_regions = {window_size: [] for window_size in args.windows}

    all_pattern_matches = []

    all_pattern_hit_seqs = []

    protein_summaries = []

    protein_count = 0

    # --------------------------------------------------------
    # Scan FASTA
    # --------------------------------------------------------

    for record in SeqIO.parse(
        args.input,
        "fasta",
    ):

        protein_count += 1

        (
            windows,
            regions,
            summary,
            pattern_matches,
            pattern_hit_seqs,
        ) = process_protein(
            record,
            args.windows,
            args.step,
            args.min_fraction,
            args.max_gap,
            args.residues,
            pattern_regex,
        )

        if summary is not None:

            protein_summaries.append(summary)

        if pattern_regex is not None:

            all_pattern_matches.extend(pattern_matches)

            all_pattern_hit_seqs.extend(pattern_hit_seqs)

        for window_size in args.windows:

            all_windows[window_size].extend(
                windows.get(
                    window_size,
                    [],
                )
            )

            all_regions[window_size].extend(
                regions.get(
                    window_size,
                    [],
                )
            )

    # --------------------------------------------------------
    # Column definitions
    # --------------------------------------------------------

    window_columns = [
        "protein_id",
        "window_size",
        "start",
        "end",
        "length",
        "selected_count",
        "fraction",
        "global_fraction",
        "enrichment",
        "longest_selected_run",
        "entropy",
        "sequence",
        "masked_sequence",
    ]

    region_columns = [
        "protein_id",
        "window_size",
        "start",
        "end",
        "length",
        "selected_count",
        "fraction",
        "global_fraction",
        "enrichment",
        "longest_selected_run",
        "longest_selected_start",
        "longest_selected_end",
        "longest_selected_sequence",
        "entropy",
        "sequence",
        "masked_sequence",
    ]

    # --------------------------------------------------------
    # Sort and write per-window / per-region files
    # --------------------------------------------------------

    window_counts = {}
    region_counts = {}

    for window_size in args.windows:

        window_rows = all_windows[window_size]

        region_rows = all_regions[window_size]

        # ----------------------------------------------------
        # Windows:
        # 1. Longest consecutive selected-residue run
        # 2. Selected-residue fraction
        # 3. Enrichment
        # ----------------------------------------------------

        window_rows.sort(
            key=lambda x: (
                x["longest_selected_run"],
                x["fraction"],
                x["enrichment"],
            ),
            reverse=True,
        )

        # ----------------------------------------------------
        # Regions:
        # 1. Enrichment
        # 2. Selected-residue fraction
        # 3. Longest consecutive selected-residue run
        # ----------------------------------------------------

        region_rows.sort(
            key=lambda x: (
                x["enrichment"],
                x["fraction"],
                x["longest_selected_run"],
            ),
            reverse=True,
        )

        window_counts[window_size] = len(window_rows)

        region_counts[window_size] = len(region_rows)

        # ----------------------------------------------------
        # Add rank to windows
        # ----------------------------------------------------

        ranked_windows = []

        for rank, row in enumerate(
            window_rows,
            start=1,
        ):

            ranked_row = row.copy()

            ranked_row["rank"] = rank

            ranked_windows.append(ranked_row)

        # ----------------------------------------------------
        # Add rank to regions
        # ----------------------------------------------------

        ranked_regions = []

        for rank, row in enumerate(
            region_rows,
            start=1,
        ):

            ranked_row = row.copy()

            ranked_row["rank"] = rank

            ranked_regions.append(ranked_row)

        # ----------------------------------------------------
        # Output files
        # ----------------------------------------------------

        window_file = output_dir / "windows" / f"{window_size}aa" / "windows.tsv"

        region_file = output_dir / "regions" / f"{window_size}aa" / "regions.tsv"

        write_tsv(
            window_file,
            ranked_windows,
            [
                "rank",
                *window_columns,
            ],
        )

        write_tsv(
            region_file,
            ranked_regions,
            [
                "rank",
                *region_columns,
            ],
        )

        print(
            f"{window_size:>3} aa | "
            f"{len(window_rows):>8} windows | "
            f"{len(region_rows):>8} regions"
        )

    # --------------------------------------------------------
    # Rank proteins by longest consecutive selected residues
    # --------------------------------------------------------

    consecutive_columns = [
        "rank",
        "protein_id",
        "length",
        "total_selected",
        "global_fraction",
        "longest_selected_run",
        "longest_selected_start",
        "longest_selected_end",
        "longest_selected_sequence",
    ]

    protein_summaries.sort(
        key=lambda x: (
            x["longest_selected_run"],
            x["global_fraction"],
            x["length"],
        ),
        reverse=True,
    )

    ranked_consecutive = []

    for rank, row in enumerate(
        protein_summaries,
        start=1,
    ):

        ranked_row = row.copy()

        ranked_row["rank"] = rank

        ranked_consecutive.append(ranked_row)

    consecutive_file = output_dir / "rankings" / "longest_consecutive_selected.tsv"

    write_tsv(
        consecutive_file,
        ranked_consecutive,
        consecutive_columns,
    )

    # --------------------------------------------------------
    # Rank proteins by global selected-residue density
    # --------------------------------------------------------

    density_columns = [
        "rank",
        "protein_id",
        "length",
        "total_selected",
        "global_fraction",
        "longest_selected_run",
        "longest_selected_start",
        "longest_selected_end",
        "longest_selected_sequence",
    ]

    ranked_density = sorted(
        protein_summaries,
        key=lambda x: (
            x["global_fraction"],
            x["longest_selected_run"],
        ),
        reverse=True,
    )

    ranked_density_output = []

    for rank, row in enumerate(
        ranked_density,
        start=1,
    ):

        ranked_row = row.copy()

        ranked_row["rank"] = rank

        ranked_density_output.append(ranked_row)

    density_file = output_dir / "rankings" / "highest_global_density.tsv"

    write_tsv(
        density_file,
        ranked_density_output,
        density_columns,
    )

    # --------------------------------------------------------
    # Pattern matches
    # --------------------------------------------------------

    pattern_match_count = 0

    if pattern_regex is not None:

        pattern_file = output_dir / "pattern_matches.tsv"

        write_tsv(
            pattern_file,
            all_pattern_matches,
            [
                "protein_id",
                "start",
                "end",
                "length",
                "sequence",
            ],
        )

        pattern_match_count = len(all_pattern_matches)

        if all_pattern_hit_seqs:

            pattern_fasta = output_dir / "pattern_hits.faa"

            with open(pattern_fasta, "w") as out:

                for hit in all_pattern_hit_seqs:

                    out.write(
                        f">{hit['protein_id']}\n"
                        f"{hit['sequence']}\n"
                    )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    summary_rows = []

    for window_size in args.windows:

        summary_rows.append(
            {
                "window_size": window_size,
                "high_fraction_windows": window_counts[window_size],
                "regions": region_counts[window_size],
            }
        )

    summary_file = output_dir / "summary.tsv"

    write_tsv(
        summary_file,
        summary_rows,
        [
            "window_size",
            "high_fraction_windows",
            "regions",
        ],
    )

    # --------------------------------------------------------
    # Log
    # --------------------------------------------------------

    end_time = datetime.now()

    elapsed = (end_time - start_time).total_seconds()

    log_file = output_dir / "run.log"

    write_log(
        log_file,
        args,
        protein_count,
        window_counts,
        region_counts,
        elapsed,
    )

    # --------------------------------------------------------
    # Final output
    # --------------------------------------------------------

    print()

    print("==================================================")

    print(f"Proteins scanned: " f"{protein_count}")

    print(f"Longest selected-residue ranking: " f"{consecutive_file}")

    print(f"Global density ranking: " f"{density_file}")

    if pattern_regex is not None:

        print(f"Pattern matches: " f"{pattern_match_count}")

        print(f"Pattern file: " f"{pattern_file}")

        if all_pattern_hit_seqs:

            print(f"Pattern hit FASTA: " f"{pattern_fasta}")

    print(f"Summary: " f"{summary_file}")

    print(f"Log: " f"{log_file}")

    print(f"Runtime: " f"{elapsed:.2f} seconds")

    print("==================================================")


if __name__ == "__main__":
    main()
