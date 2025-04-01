import csv
from astropy.io.votable import parse
from collections import defaultdict

# Load the VOTable file
votable = parse("Variable_objects/Welch_Stetson_Variables/matched_variables_ws.vot")  # Change when necessary

# Access the first table in the VOTable
table = votable.get_first_table().to_table()

# Extract relevant columns
stellar_types = table["OTYPE_S"]  # Stellar types (OTYPE_S)
stetson_values = table["Welch_Stetson"]  # Stetson values

# Dictionary to store sum of Stetson values and count per OTYPE_S
stellar_stats = defaultdict(lambda: {"sum": 0, "count": 0})

# Accumulate sums and counts
for stellar_type, stetson in zip(stellar_types, stetson_values):
    stellar_stats[stellar_type]["sum"] += stetson
    stellar_stats[stellar_type]["count"] += 1

# Save results to CSV file
csv_filename = "otype_av.csv"
with open(csv_filename, mode="w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["OTYPE_S", "Count", "Avg_Index"])  # CSV header
    
    for stellar_type, stats in stellar_stats.items():
        avg_stetson = stats["sum"] / stats["count"]
        writer.writerow([stellar_type, stats["count"], round(avg_stetson, 4)])

print(f"CSV file '{csv_filename}' saved successfully!")

