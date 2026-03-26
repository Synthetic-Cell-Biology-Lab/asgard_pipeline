import sys
import os
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# -----------------------------
# Input
# -----------------------------
fasta = sys.argv[1]
output = sys.argv[2]
split_dir = sys.argv[3]


records = list(SeqIO.parse(fasta, "fasta"))

lengths = [len(r.seq) for r in records]
df = pd.DataFrame({
    "Length": lengths,
    "Record": records
})


df["quantile_bin"] = pd.qcut(df["Length"], q=4, labels = ["short", "mid_short", "mid_long", "long"])

os.makedirs(split_dir, exist_ok=True)

for q in df["quantile_bin"].unique():
    subset = df[df["quantile_bin"] == q]

    out_path = os.path.join(split_dir, f"{q}.fasta")

    SeqIO.write(subset["Record"].tolist(), out_path, "fasta")

    print(f"{q}: {len(subset)} sequences → {out_path}")

# -----------------------------
# Style (clean + aesthetic)
# -----------------------------
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams["figure.dpi"] = 120

# -----------------------------
# Create layout
# -----------------------------
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.25)

# -----------------------------
# 1. Histogram
# -----------------------------
ax1 = fig.add_subplot(gs[0, 0])
sns.histplot(df["Length"], bins=40, kde=True, ax=ax1)
ax1.set_title("Histogram + KDE")
ax1.set_xlabel("Sequence Length")
ax1.set_ylabel("Frequency")

quantiles = df["Length"].quantile([0.25, 0.5, 0.75])

for q in quantiles:
    ax1.axvline(q, linestyle="--", alpha=0.7)


# -----------------------------
# 2. Violin plot
# -----------------------------
ax2 = fig.add_subplot(gs[0, 1])
sns.violinplot(y=df["Length"], inner="quartile", ax=ax2)
ax2.set_title("Violin Plot")
ax2.set_ylabel("Sequence Length")

# -----------------------------
# 3. Bar plot (binned counts)
# -----------------------------
ax3 = fig.add_subplot(gs[1, 0])
bins = np.linspace(df["Length"].min(), df["Length"].max(), 15)
df["bin"] = pd.cut(df["Length"], bins)

bin_counts = df["bin"].value_counts().sort_index()
bin_labels = [f"{int(b.left)}-{int(b.right)}" for b in bin_counts.index]

ax3.bar(bin_labels, bin_counts.values)
ax3.set_title("Binned Counts")
ax3.set_xlabel("Length Range")
ax3.set_ylabel("Count")
ax3.tick_params(axis='x', rotation=45)

# -----------------------------
# 4. Beeswarm (strip plot)
# -----------------------------
ax4 = fig.add_subplot(gs[1, 1])
sns.swarmplot(y=df["Length"], size=3, ax=ax4)
ax4.set_title("Beeswarm Plot")
ax4.set_ylabel("Sequence Length")

# -----------------------------
# Global title
# -----------------------------
fig.suptitle("Sequence Length Distribution Overview", fontsize=18)

# -----------------------------
# Save
# -----------------------------
plt.savefig(output, bbox_inches="tight")
plt.close()


df[["Length", "quantile_bin"]].to_csv("length_bins.csv", index=False)