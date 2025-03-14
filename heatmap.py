import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("count_object_type.csv")

# Stellar taxonomy (same as csv file headers) + total number of objects
object_totals = {
    "YSOs": 2918,
    "MS stars": 9,
    "Evolved stars": 332,
    "Interacting binaries": 1157,
    "Variability" : 5,
    "General spectral properties" : 491,
    "Star": 1468
}

df["Fraction"] = df.apply(
    lambda row: row["Count"] / object_totals[row["Object_type"]],
    axis=1
)

heatmap_data = df.pivot(index="Object_type", columns="Index", values="Fraction")

plt.figure(figsize=(11, 7))
ax = sns.heatmap(
    heatmap_data, 
    cmap="viridis", 
    annot=True, 
    fmt=".2f", 
    cbar_kws={"label": "Ratio detected"}
)

new_labels = [f"{label} ({object_totals[label]})" for label in heatmap_data.index] # Include total number of objects in y ticks
ax.set_yticklabels(new_labels, rotation=0)

plt.title("Fraction of Each Object Type Flagged by Each Variability Index", fontweight="bold", size=12)
plt.xlabel("Variability Index", fontweight="bold", size=12)
plt.ylabel("Object Type (Total nº of objects)", fontweight="bold", size=12)
plt.xticks(size=10, color='black')
plt.tight_layout()
plt.savefig('Figures/heatmap.png')
plt.show()

