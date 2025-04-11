import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("object_type_count.csv")

object_totals = { 
    "YSO": 2920,
    "MS star": 9,
    "Evolved stars": 332,
    "Interacting binaries": 1157,
    "Substellar objects": 3,
    "Environment properties": 34,
    "Star": 1468
} # Total number of objects in the SIMBAD database for the reigon

df["Fraction"] = df.apply(
    lambda row: row["Count"] / object_totals[row["Object_type"]],
    axis=1
)  # Calculates ratio detected 

heatmap_data = df.pivot(index="Object_type", columns="Index", values="Fraction")

plt.figure(figsize=(12, 8))
ax = sns.heatmap(
    heatmap_data, 
    cmap="viridis", 
    annot=True, 
    fmt=".2f", 
    cbar_kws={"label": "Ratio detected"}
)
ax.collections[0].colorbar.set_label("Ratio detected", size=13, weight="bold")

new_labels = [f"{label} ({object_totals[label]})" for label in heatmap_data.index]
ax.set_yticklabels(new_labels, rotation=0)

# Set plot
plt.title("Fraction of Each Object Type Flagged by Each Variability Index", fontweight="bold", size=13)
plt.xlabel("Variability Index", fontweight="bold", size=15)
plt.ylabel("Object Type (Total nº of objects)", fontweight="bold", size=15)
plt.xticks(size=12, color='black')
plt.yticks(size=12, color='black')
plt.tight_layout()
plt.savefig('Figures/heatmap.png')
plt.show()

