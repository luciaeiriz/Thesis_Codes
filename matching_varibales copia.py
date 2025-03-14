from astropy.io.votable import parse, from_table, writeto
import os
from astropy.table import Table

# Folder containing the VOTables
folder = "Variable_objects/SA04_Variables" #change when necessary

# List to store tables
tables = []

# Read all VOTables in the folder
for filename in os.listdir(folder):
    if filename.endswith(".vot"):  #change
        votable = parse(os.path.join(folder, filename))
        table = votable.get_first_table().to_table()
        tables.append(table)

# Extract unique stellar object identifiers (RA, DEC) from each table
stellar_counts = {}

for table in tables:
    for row in table:
        ra_dec = (row[0], row[1])  # Assuming RA is first column, DEC is second column
        stellar_counts[ra_dec] = stellar_counts.get(ra_dec, 0) + 1

# Filter objects that appear in more than three tables
selected_ids = {ra_dec for ra_dec, count in stellar_counts.items() if count >= 3}

# Create a new VOTable with the matched objects
matched_table = Table(names=tables[0].colnames, dtype=[tables[0][col].dtype for col in tables[0].colnames])

# Set to keep track of added stars
added_stars = set()

# Fill the new table with rows matching the selected RA, DEC pairs, avoiding duplicates
for table in tables:
    for row in table:
        ra_dec = (row[0], row[1])
        if ra_dec in selected_ids and ra_dec not in added_stars:
            matched_table.add_row(row)
            added_stars.add(ra_dec)


# Create a new VOTable and save it
index = "sa04" #change when necessary
output_path = os.path.join(folder, "variables_"+index+".vot")
new_votable = from_table(matched_table)
writeto(new_votable, output_path)
print(f"New VOTable 'matched_stars.vot' created successfully in {folder}")
