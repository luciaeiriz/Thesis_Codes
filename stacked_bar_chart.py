import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read in your CSV file
df = pd.read_csv("count_object_type.csv")

# Ensure that the 'Object_type' column has no extra whitespace
df['Object_type'] = df['Object_type'].str.strip()

# Define object totals
object_totals = {
    "YSOs": 2918,
    "MS stars": 9,
    "Evolved stars": 332,
    "Interacting binaries": 1157,
    "Variability" : 5,
    "General spectral properties" : 491,
    "Star": 1468
}

# Pivot the DataFrame to organize data for a stacked/grouped bar chart
bar_data = df.pivot(index="Index", columns="Object_type", values="Count")

# Create figure and axis
fig, ax = plt.subplots(figsize=(12, 7))

# Define a colormap based on "viridis"
colors = sns.color_palette("viridis", n_colors=len(bar_data.columns))

# ---- Option 1: Stacked Bar Chart ----
# bar_data.plot(kind="bar", stacked=True, color=colors, edgecolor="black", alpha=0.85, ax=ax)

# ---- Option 2: Grouped Bar Chart ----
bar_data.plot(kind="bar", stacked=False, color=colors, edgecolor="black", alpha=0.85, width=0.8, position=1, ax=ax)

# Modify legend labels to include total counts
new_labels = [f"{label} ({object_totals[label]})" for label in bar_data.columns]
ax.legend(new_labels, title="Object Type (Total Count)", fontsize=10, bbox_to_anchor=(1.05, 1), loc="upper left")

# Formatting
plt.title("Count of Each Object Type Flagged by Variability Indices", fontsize=14, fontweight="bold")
plt.xlabel("Variability Index", fontsize=12, fontweight="bold")
plt.ylabel("Count", fontsize=12, fontweight="bold")
plt.xticks(rotation=45, ha="right", fontsize=10)

plt.tight_layout()
plt.savefig('Figures/bar_chart.png')
plt.show()

