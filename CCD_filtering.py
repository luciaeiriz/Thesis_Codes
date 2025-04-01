from astropy.table import Table
import numpy as np 

# Load your VOTable
folder_path = "Variable_objects/Stetson_Variables/"
table = Table.read(folder_path + "variables_stetson.vot", format="votable")

# Define filtering conditions
ra_values_to_remove = [(343.55, 62.66, 63.36), (344.45, 62.68, 63.36), (344.44, 62.68, 63.36)]
tolerance = 1e-2  # RA precision tolerance

# Create mask to exclude specified RA and DEC ranges
mask = np.ones(len(table), dtype=bool)  # Start with all True (keep all rows)

for ra, dec_min, dec_max in ra_values_to_remove:
    match_ra = np.isclose(table['RA'], ra, atol=tolerance)
    match_dec = (table['DEC'] >= dec_min) & (table['DEC'] <= dec_max)
    mask &= ~(match_ra & match_dec)  # Exclude matching rows

# Apply the mask
filtered_table = table[mask]

# Save the filtered table
filtered_table.write(folder_path + "filtered_CCD_variables.vot", format="votable", overwrite=True)

print("Filtering complete! Saved as 'filtered_CCD_variables.vot'.")


