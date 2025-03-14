from astropy.io.votable import parse
from collections import Counter

# Load the VOTable file
votable = parse("Variable_objects/SA04_Variables/simbad_match_sa04.vot")  # change when necessary

# Access the first table in the VOTable
table = votable.get_first_table().to_table()

# Extract the first column (stellar type)
stellar_types = table[table.colnames[0]]  # Assuming the first column has stellar types

# Count occurrences of each stellar type
stellar_counts = Counter(stellar_types)

# Print results
for stellar_type, count in stellar_counts.items():
    print(f"{stellar_type} {count}")
