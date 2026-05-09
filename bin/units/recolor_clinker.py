#!/usr/bin/env python3

import argparse
import json
import re
import uuid
from collections import defaultdict

import pandas as pd


############################################################
# Extract clinker data from HTML
############################################################

def extract_data_from_html(html_text):

    pattern = re.compile(
        r"const\s+data\s*=\s*(\{.*?\})\s*;function",
        re.DOTALL
    )

    match = pattern.search(html_text)

    if not match:
        raise ValueError(
            "Could not find embedded clinker data block"
        )

    json_text = match.group(1)

    return json.loads(json_text)


############################################################
# Inject modified data back into HTML
############################################################

def inject_data_into_html(html_text, session):

    new_json = json.dumps(session)

    pattern = re.compile(
        r"(const\s+data\s*=\s*)(\{.*?\})(\s*;function)",
        re.DOTALL
    )

    html_text = pattern.sub(
        rf"\1{new_json}\3",
        html_text
    )

    return html_text


############################################################
# Load ontology color map
############################################################

def load_color_map(color_csv):

    df = pd.read_csv(color_csv)

    df["locus_tag"] = (
        df["locus_tag"]
        .astype(str)
        .str.strip()
    )

    df["function"] = (
        df["function"]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    df["color"] = (
        df["color"]
        .astype(str)
        .str.strip()
    )

    return df


############################################################
# Extract UID <-> locus mappings
############################################################

def extract_gene_uid_map(session):

    uid_to_locus = {}

    clusters = session.get("clusters", {})

    if isinstance(clusters, dict):
        cluster_iter = clusters.values()
    else:
        cluster_iter = clusters

    for cluster in cluster_iter:

        if not isinstance(cluster, dict):
            continue

        for locus in cluster.get("loci", []):

            if not isinstance(locus, dict):
                continue

            for gene in locus.get("genes", []):

                if not isinstance(gene, dict):
                    continue

                uid = gene.get("uid")

                ################################################
                # Prefer real locus_tag
                ################################################

                locus_tag = (
                    gene.get("names", {})
                    .get("locus_tag")
                )

                ################################################
                # Fallback
                ################################################

                if not locus_tag:
                    locus_tag = gene.get("label")

                if uid and locus_tag:

                    uid_to_locus[uid] = (
                        str(locus_tag).strip()
                    )

    return uid_to_locus


############################################################
# Build custom ontology groups
############################################################

def build_custom_groups(
    uid_to_locus,
    color_df
):

    ########################################################
    # Build ontology maps
    ########################################################

    locus_to_function = dict(
        zip(
            color_df["locus_tag"],
            color_df["function"]
        )
    )

    function_to_color = dict(
        zip(
            color_df["function"],
            color_df["color"]
        )
    )

    ########################################################
    # Group gene UIDs by ontology
    ########################################################

    function_to_genes = defaultdict(list)

    for uid, locus in uid_to_locus.items():

        function = locus_to_function.get(locus)

        if not function:
            function = "unknown"

        function_to_genes[function].append(uid)

    ########################################################
    # Build clinker group objects
    ########################################################

    new_groups = []

    for function, gene_uids in function_to_genes.items():

        color = function_to_color.get(
            function,
            "#B0B0B0"
        )

        group = {
            "uid": str(uuid.uuid4()),
            "label": function,
            "genes": gene_uids,
            "hidden": False,
            "colour": color,
            "groupColour": color,
            "color": color
        }

        new_groups.append(group)

    return new_groups


############################################################
# MAIN
############################################################

def main():

    parser = argparse.ArgumentParser(
        description="Build custom clinker ontology groups"
    )

    parser.add_argument(
        "--html",
        required=True,
        help="Original clinker HTML"
    )

    parser.add_argument(
        "--color_map",
        required=True,
        help="CSV with locus_tag,function,color"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output recolored HTML"
    )

    args = parser.parse_args()

    ########################################################
    # Load HTML
    ########################################################

    with open(args.html) as f:
        html_text = f.read()

    ########################################################
    # Extract clinker session
    ########################################################

    session = extract_data_from_html(html_text)

    ########################################################
    # Load ontology
    ########################################################

    color_df = load_color_map(args.color_map)

    ########################################################
    # Extract mappings
    ########################################################

    uid_to_locus = extract_gene_uid_map(session)

    print(
        f"Extracted {len(uid_to_locus)} genes"
    )

    ########################################################
    # Build ontology groups
    ########################################################

    new_groups = build_custom_groups(
        uid_to_locus,
        color_df
    )

    print(
        f"Built {len(new_groups)} ontology groups"
    )

    ########################################################
    # Replace clinker groups
    ########################################################

    session["groups"] = new_groups

    ########################################################
    # Inject modified data back into HTML
    ########################################################

    modified_html = inject_data_into_html(
        html_text,
        session
    )

    ########################################################
    # Save
    ########################################################

    with open(args.output, "w") as f:
        f.write(modified_html)

    print("\nSaved:")
    print(args.output)


if __name__ == "__main__":
    main()

