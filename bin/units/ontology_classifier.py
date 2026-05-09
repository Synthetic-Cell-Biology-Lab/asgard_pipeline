#!/usr/bin/env python3

import re
import argparse
import pandas as pd


############################################################
# Functional classifier rules
############################################################

CLASSIFIER_RULES = {

    ########################################################
    # Cytoskeleton
    ########################################################

    "cytoskeleton": [
        r"tubulin",
        r"ftsz",
        r"gelsolin",
        r"actin",
        r"cytoskeleton",
    ],

    ########################################################
    # Ribosome / Translation
    ########################################################

    "ribosomal": [
        r"ribosomal",
        r"\buL\d",
        r"\buS\d",
        r"\beL\d",
        r"\beS\d",
        r"ribosome",
        r"large ribosomal subunit",
        r"small ribosomal subunit",
    ],

    "translation": [
        r"translation",
        r"elongation factor",
        r"initiation factor",
        r"release factor",
        r"aminoacyl",
        r"synthetase",
        r"tRNA",
        r"rRNA",
        r"nusg",
        r"spt5",
        r"kow",
        r"ngn domain",
        r"if6",
    ],

    ########################################################
    # DNA / RNA
    ########################################################

    "helicase": [
        r"helicase",
        r"dead",
        r"deah",
        r"uvrb",
        r"radA",
        r"recA",
    ],

    "dna_binding": [
        r"dna-binding",
        r"winged helix",
        r"helix-turn-helix",
        r"\bhth\b",
        r"marr",
        r"cro",
        r"transcription regulator",
    ],

    "rna_processing": [
        r"ribonuclease",
        r"rna",
        r"pseudouridine",
        r"methyltransferase",
        r"crm domain",
        r"splice",
        r"alba",
        r"rpp21",
        r"rpr2",
        r"snm1",
    ],

    ########################################################
    # ATPases / GTPases
    ########################################################

    "atpase_gtpase": [
        r"gtpase",
        r"atpase",
        r"p-loop",
        r"\baaa\b",
        r"obg",
        r"cobw",
        r"hypb",
        r"ureg",
        r"ftsY",
        r"kaiC",
        r"small gtpase",
        r"abc transporter",
    ],

    ########################################################
    # Transport / Membrane
    ########################################################

    "transport_membrane": [
        r"transporter",
        r"transmembrane",
        r"membrane",
        r"permease",
        r"channel",
        r"sec",
        r"export",
        r"import",
        r"secretion",
        r"type ii secretion",
        r"type iv secretion",
        r"mfs transporter",
    ],

    ########################################################
    # Redox / Fe-S / Electron transfer
    ########################################################

    "redox": [
        r"ferredoxin",
        r"oxidoreductase",
        r"dehydrogenase",
        r"electron transfer",
        r"redox",
        r"thioredoxin",
        r"nitroreductase",
        r"flavoprotein",
        r"radical sam",
        r"iron-sulphur",
        r"fe-s",
        r"nif",
    ],

    ########################################################
    # Metal / Cofactor systems
    ########################################################

    "metal_cofactor": [
        r"metal",
        r"cobalamin",
        r"vitamin b12",
        r"nickel",
        r"copper",
        r"iron",
        r"fe-s",
        r"nfu",
        r"hyp",
        r"cob",
        r"arsenate",
    ],

    ########################################################
    # Metabolism
    ########################################################

    "metabolism": [
        r"kinase",
        r"epimerase",
        r"hydrolase",
        r"transferase",
        r"deaminase",
        r"synthase",
        r"carboxylase",
        r"biosynthesis",
        r"phosphatase",
        r"amidase",
        r"reductase",
        r"ligase",
    ],

    ########################################################
    # Chaperones / Folding
    ########################################################

    "protein_folding": [
        r"chaperone",
        r"dnak",
        r"hsp70",
        r"prefoldin",
        r"folding",
    ],

    ########################################################
    # Stress / Repair
    ########################################################

    "stress_repair": [
        r"dna repair",
        r"repair",
        r"stress",
        r"damage",
        r"dinb",
    ],

    ########################################################
    # Regulatory / Signaling
    ########################################################

    "regulatory": [
        r"sensor",
        r"signaling",
        r"regulator",
        r"response regulator",
        r"transcription factor",
    ],

    ########################################################
    # Unknown conserved
    ########################################################

    "unknown": [
        r"duf",
        r"uncharacterised",
        r"hypothetical protein",
        r"domain of unknown function",
    ],
}


############################################################
# Compile regex patterns
############################################################

COMPILED_RULES = {
    category: [
        re.compile(pattern, re.IGNORECASE)
        for pattern in patterns
    ]
    for category, patterns in CLASSIFIER_RULES.items()
}


############################################################
# Category priority
############################################################

CATEGORY_PRIORITY = [
    "cytoskeleton",
    "ribosomal",
    "translation",
    "helicase",
    "dna_binding",
    "rna_processing",
    "atpase_gtpase",
    "transport_membrane",
    "redox",
    "metal_cofactor",
    "protein_folding",
    "stress_repair",
    "regulatory",
    "metabolism",
    "unknown",
]


############################################################
# Classify annotation
############################################################

def classify_annotation(annotation):

    if pd.isna(annotation):
        return "unknown"

    annotation = str(annotation)

    scores = {}

    for category, patterns in COMPILED_RULES.items():

        score = 0

        for pattern in patterns:

            if pattern.search(annotation):
                score += 1

        scores[category] = score

    ########################################################
    # Highest score wins
    ########################################################

    best_score = max(scores.values())

    if best_score == 0:
        return "other"

    ########################################################
    # Resolve ties using priority
    ########################################################

    tied = [
        category
        for category, score in scores.items()
        if score == best_score
    ]

    for category in CATEGORY_PRIORITY:
        if category in tied:
            return category

    return "other"


############################################################
# Apply classification
############################################################

def classify_dataframe(df, annotation_column):

    df["functional_bin"] = df[annotation_column].apply(
        classify_annotation
    )

    return df


############################################################
# MAIN
############################################################

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--column", required=True)

    ########################################################
    # Outputs
    ########################################################

    parser.add_argument(
        "--classified_output",
        required=True,
        help="Full classified dataframe output"
    )


    args = parser.parse_args()

    ########################################################
    # Load dataframe
    ########################################################

    df = pd.read_csv(args.input)

    ########################################################
    # Classify
    ########################################################

    df = classify_dataframe(
        df,
        annotation_column=args.column
    )

    ########################################################
    # Save full dataframe
    ########################################################

    df.to_csv(
        args.classified_output,
        index=False
    )

    ########################################################
    # Summary
    ########################################################

    print("\nFunctional bin counts:\n")
    print(df["functional_bin"].value_counts())

    print(f"\nSaved classified dataframe:")
    print(args.classified_output)


if __name__ == "__main__":
    main()