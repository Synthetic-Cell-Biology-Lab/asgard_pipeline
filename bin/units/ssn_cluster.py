#!/usr/bin/env python3
"""
SSN Step 6: Cluster assignment + per-cluster FASTA files.

For each bitscore level, computes connected components from the edge list,
writes a clusters CSV (one column per bitscore), and splits the NR FASTA
into per-cluster FASTA files.

Usage:
    python ssn_cluster.py \
        --nodes     <nodes.tsv>           \
        --fasta     <nr.fasta>            \
        --edges     <bs50.edges.tsv> <bs100.edges.tsv> ... \
        --bitscores 50,100,150            \
        --out-csv   <clusters.csv>        \
        --out-dir   <cluster_fastas/>
"""

import argparse
import csv
import logging
import os
import sys
from collections import defaultdict

# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------

logging.basicConfig(
    format="[%(asctime)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
    stream=sys.stdout,
)
log = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Union-Find
# -----------------------------------------------------------------------------

class UnionFind:
    def __init__(self):
        self._parent: dict[str, str] = {}

    def find(self, x: str) -> str:
        p = self._parent
        while p.setdefault(x, x) != x:
            p[x] = p[p[x]]   # path compression
            x = p[x]
        return x

    def union(self, a: str, b: str) -> None:
        self._parent[self.find(a)] = self.find(b)

    def components(self) -> dict[str, int]:
        """Return {node: cluster_label} with labels assigned in sorted order."""
        root_to_label: dict[str, int] = {}
        counter = 1
        mapping: dict[str, int] = {}
        for node in sorted(self._parent):
            root = self.find(node)
            if root not in root_to_label:
                root_to_label[root] = counter
                counter += 1
            mapping[node] = root_to_label[root]
        return mapping


# -----------------------------------------------------------------------------
# I/O helpers
# -----------------------------------------------------------------------------

def load_node_ids(nodes_tsv: str) -> list[str]:
    with open(nodes_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [row["id"] for row in reader]


def load_fasta(fasta_path: str) -> dict[str, tuple[str, str]]:
    """Return {seq_id: (header_line, sequence)}."""
    sequences: dict[str, tuple[str, str]] = {}
    current_id: str | None = None
    current_header: str = ""
    current_seq: list[str] = []

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = (current_header, "".join(current_seq))
                current_header = line
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    if current_id is not None:
        sequences[current_id] = (current_header, "".join(current_seq))

    return sequences


def connected_components(edge_file: str) -> dict[str, int]:
    uf = UnionFind()
    with open(edge_file) as fh:
        for line in fh:
            a, b, *_ = line.strip().split("\t")
            uf.union(a, b)
    return uf.components()


def write_fasta(path: str, members: list[str], sequences: dict[str, tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for node in members:
            if node not in sequences:
                log.warning("Node %s not found in FASTA — skipping.", node)
                continue
            header, seq = sequences[node]
            fh.write(f"{header}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--nodes",     required=True, help="nodes.tsv from ssn_filter")
    p.add_argument("--fasta",     required=True, help="NR FASTA from ssn_cdhit")
    p.add_argument("--edges",     required=True, nargs="+", help="one edge TSV per bitscore level")
    p.add_argument("--bitscores", required=True, help="comma-separated bitscore values matching --edges order")
    p.add_argument("--out-csv",   required=True, help="output clusters CSV")
    p.add_argument("--out-dir",   required=True, help="output directory for per-cluster FASTAs")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    bitscores = [b.strip() for b in args.bitscores.split(",")]

    if len(bitscores) != len(args.edges):
        log.error(
            "Number of bitscore values (%d) does not match number of edge files (%d).",
            len(bitscores), len(args.edges),
        )
        sys.exit(1)

    # ── Load inputs ───────────────────────────────────────────────────────────

    log.info("Loading node IDs from %s", args.nodes)
    all_nodes = load_node_ids(args.nodes)
    log.info("  %d nodes", len(all_nodes))

    log.info("Parsing FASTA from %s", args.fasta)
    sequences = load_fasta(args.fasta)
    log.info("  %d sequences", len(sequences))

    # ── Connected components per bitscore ─────────────────────────────────────

    bs_mappings: dict[str, dict[str, int]] = {}
    for bs, ef in zip(bitscores, args.edges):
        log.info("Computing connected components for bitscore=%s from %s", bs, ef)
        bs_mappings[bs] = connected_components(ef)
        n_clusters = len(set(bs_mappings[bs].values()))
        n_connected = len(bs_mappings[bs])
        log.info("  %d clusters, %d connected nodes", n_clusters, n_connected)

    # ── Write clusters CSV ────────────────────────────────────────────────────

    os.makedirs(os.path.dirname(args.out_csv) or ".", exist_ok=True)
    log.info("Writing clusters CSV -> %s", args.out_csv)

    fieldnames = ["id"] + [f"cluster_bs{bs}" for bs in bitscores]
    with open(args.out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for node in all_nodes:
            row: dict[str, str | int] = {"id": node}
            for bs in bitscores:
                row[f"cluster_bs{bs}"] = bs_mappings[bs].get(node, 0)
            writer.writerow(row)

    # ── Write per-cluster FASTAs ──────────────────────────────────────────────

    log.info("Writing per-cluster FASTAs -> %s", args.out_dir)

    for bs in bitscores:
        bs_dir = os.path.join(args.out_dir, f"bs{bs}")
        os.makedirs(bs_dir, exist_ok=True)

        clusters: dict[int, list[str]] = defaultdict(list)
        for node in all_nodes:
            label = bs_mappings[bs].get(node, 0)
            clusters[label].append(node)

        for label, members in sorted(clusters.items()):
            fname = "singletons.fasta" if label == 0 else f"cluster_{label}.fasta"
            write_fasta(os.path.join(bs_dir, fname), members, sequences)

        n_real = sum(1 for l in clusters if l != 0)
        n_singletons = len(clusters.get(0, []))
        log.info("  bs%s: %d clusters, %d singletons", bs, n_real, n_singletons)

    log.info("Cluster step completed successfully.")


if __name__ == "__main__":
    main()