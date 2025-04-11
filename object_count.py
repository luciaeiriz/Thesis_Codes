from astropy.io.votable import parse
from collections import Counter

# Load the VOTable file
votable = parse("Variable_objects/SA04_Variables/matched_variables_sa04.vot")  # change when necessary

# Access the first table in the VOTable
table = votable.get_first_table().to_table()

# Extract relevant columns
stellar_types = table[table.colnames[0]]  # First column has stellar types
ra_values = table[table.colnames[4]]  # Second column is RA
dec_values = table[table.colnames[5]]  # Third column is DEC

# Count occurrences of each stellar type
stellar_counts = Counter(stellar_types)

# Print results with corresponding RA_d and DEC_d
for stellar_type, count in stellar_counts.items():
    print(f"{stellar_type} {count}")
    # print("RA, DEC:")
    # for i in range(len(stellar_types)):
    #     if stellar_types[i] == stellar_type:
    #         print(f"{ra_values[i]}, {dec_values[i]}")
    # print("-")

