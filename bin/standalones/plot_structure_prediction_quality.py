import os
import numpy as np
import matplotlib.pyplot as plt
import pymol

# ---- SETTINGS ----
Protein_folder = "/home/anirudh/asgard_pipeline/database/structures/cetz/"   # folder containing PDB files
pdb_folder = os.path.join(Protein_folder, "pdbs")
n_bins = 100

all_plddt = []
all_ss = []

with pymol2.PyMOL() as pymol:
    cmd = pymol.cmd

    for file in os.listdir(pdb_folder):
        if not file.endswith(".pdb"):
            continue

        path = os.path.join(pdb_folder, file)
        obj_name = "prot"

        try:
            cmd.load(path, obj_name)
            cmd.remove("not name CA")   # only CA atoms

            # assign secondary structure
            cmd.dss(obj_name)

            model = cmd.get_model(obj_name)

            plddt = []
            ss = []

            for atom in model.atom:
                plddt.append(atom.b)  # pLDDT
                s = atom.ss

                # map SS to categories
                if s == 'H':
                    ss.append(0)  # helix
                elif s == 'S':
                    ss.append(1)  # sheet
                else:
                    ss.append(2)  # coil

            cmd.delete(obj_name)

        except Exception as e:
            print(f"Skipping {file}: {e}")
            cmd.delete("all")
            continue

        if len(plddt) < 10:
            continue

        plddt = np.array(plddt)
        ss = np.array(ss)

        # ---- BINNING ----
        positions = np.linspace(0, 1, len(plddt))
        bin_edges = np.linspace(0, 1, n_bins + 1)

        plddt_bins = np.zeros(n_bins)
        counts = np.zeros(n_bins)
        ss_bins = np.zeros((n_bins, 3))

        for pos, p, s in zip(positions, plddt, ss):
            idx = np.searchsorted(bin_edges, pos) - 1
            idx = min(max(idx, 0), n_bins - 1)

            plddt_bins[idx] += p
            counts[idx] += 1
            ss_bins[idx, s] += 1

        valid = counts > 0
        plddt_bins[valid] /= counts[valid]
        ss_bins[valid] /= counts[valid][:, None]

        all_plddt.append(plddt_bins)
        all_ss.append(ss_bins)

# ---- AGGREGATE ----
all_plddt = np.array(all_plddt)
all_ss = np.array(all_ss)

mean_plddt = np.nanmean(all_plddt, axis=0)
mean_ss = np.nanmean(all_ss, axis=0)

x = np.linspace(0, 1, n_bins)

# ---- PLOT ----
fig, ax1 = plt.subplots(figsize=(9, 5))

ax1.plot(x, mean_plddt, label="Mean pLDDT")
ax1.set_xlabel("Normalized Position")
ax1.set_ylabel("pLDDT")

ax2 = ax1.twinx()
ax2.plot(x, mean_ss[:, 0], linestyle="--", label="Helix")
ax2.plot(x, mean_ss[:, 1], linestyle="--", label="Sheet")
ax2.plot(x, mean_ss[:, 2], linestyle="--", label="Coil")
ax2.set_ylabel("SS Fraction")

lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2)

plt.title("pLDDT vs Secondary Structure (PyMOL)")
plt.tight_layout()

img_out = os.path.join(Protein_folder, "pLDDT_profile.png")
plt.savefig(img_out, dpi=300)
plt.show()