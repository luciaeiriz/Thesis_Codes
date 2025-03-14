from astropy.table import Table
import numpy as np 

# Load your VOTable
table = Table.read("Variable_objects/Chi_Squared_Variables/variables_chisquared.vot", format="votable")

# Define RA and DEC filtering conditions with a tolerance
ra_values_to_remove = [343.36, 343.55]
dec_min, dec_max = 62.27, 63.38
tolerance = 1e-2  # Adjust based on decimal precision

# Create a mask using np.isclose for RA and range check for DEC
mask = ~((np.any([np.isclose(table['RA'], ra, atol=tolerance) for ra in ra_values_to_remove], axis=0)) &
         (table['DEC'] >= dec_min) & (table['DEC'] <= dec_max))

# Apply the mask
filtered_table = table[mask]

# Save the filtered table
filtered_table.write("filtered_file.vot", format="votable", overwrite=True)

print("Filtering complete! Saved as 'filtered_file.vot'.")

