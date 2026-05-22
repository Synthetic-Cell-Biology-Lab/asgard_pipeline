# script to give the secondary structure information given the tertiary structure of proteins in cif format

from pathlib import Path
from Bio.PDB import MMCIFParser
from Bio.PDB.DSSP import DSSP

# Folder containing CIF files
folder_path = "/home/anirudh/asgard_pipeline/database/structures/comparisons/Target"

# Parse all CIF files recursively
cif_files = list(Path(folder_path).rglob("*.cif"))

parser = MMCIFParser(QUIET=True)

for cif_file in cif_files:
    print(f"\nProcessing: {cif_file}")

    try:
        structure = parser.get_structure(cif_file.stem, str(cif_file))
        model = structure[0]

        # Run DSSP
        dssp = DSSP(model, str(cif_file), file_type="MMCIF")

        # Build linear secondary structure string
        linear_ss = []

        for key in dssp.keys():
            ss = dssp[key][2]

            # Simplify DSSP classes
            if ss in ["H", "G", "I"]:
                linear_ss.append("H")   # Helix
            elif ss in ["E", "B"]:
                linear_ss.append("E")   # Beta sheet
            else:
                linear_ss.append("C")   # Coil/Loop

        linear_ss = "".join(linear_ss)

        print("Secondary structure:")
        print(linear_ss)

    except Exception as e:
        print(f"Failed: {e}")