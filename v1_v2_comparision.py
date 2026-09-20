import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Get CSV paths
csv1_path = input("Enter path to first CSV: ").strip()
csv2_path = input("Enter path to second CSV: ").strip()

# Read CSV files
df1 = pd.read_csv(csv1_path)
df2 = pd.read_csv(csv2_path)

# Check class column
if "class" not in df1.columns or "class" not in df2.columns:
    raise ValueError("Both CSV files must contain a 'class' column.")

# Clean class names
class1 = df1["class"].dropna().astype(str).str.strip()
class2 = df2["class"].dropna().astype(str).str.strip()

class1 = class1[class1 != ""]
class2 = class2[class2 != ""]

# Count organisms in each class
counts1 = class1.value_counts()
counts2 = class2.value_counts()

# Combine all classes from both datasets
all_classes = sorted(set(counts1.index) | set(counts2.index))

# Get counts, using 0 when a class is absent
values1 = [counts1.get(c, 0) for c in all_classes]
values2 = [counts2.get(c, 0) for c in all_classes]

# X positions
x = np.arange(len(all_classes))
width = 0.4

# Create graph
fig, ax = plt.subplots(figsize=(14, 8))

bars1 = ax.bar(
    x - width / 2,
    values1,
    width,
    label="Version 1"
)

bars2 = ax.bar(
    x + width / 2,
    values2,
    width,
    label="Version 2"
)

# Add count labels
ax.bar_label(
    bars1,
    labels=[str(v) for v in values1],
    padding=3,
    fontsize=9
)

ax.bar_label(
    bars2,
    labels=[str(v) for v in values2],
    padding=3,
    fontsize=9
)

# Labels and title
ax.set_xlabel("Taxonomic Class")
ax.set_ylabel("Number of Organisms")
ax.set_title("Comparison of Organism Counts by Class")

ax.set_xticks(x)
ax.set_xticklabels(all_classes, rotation=45, ha="right")

ax.legend()

plt.tight_layout()

# Save
plt.savefig(
    "class_comparison.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("\nDataset 1:")
print(counts1)

print("\nDataset 2:")
print(counts2)

print("\nGraph saved as: class_comparison.png")