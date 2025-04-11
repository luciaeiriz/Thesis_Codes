from astropy.table import Table
import numpy as np 

#Folder path
folder_path = "Variable_objects/Stetson_Variables/"
table = Table.read(folder_path + "variables_stetson.vot", format="votable")

#Coordinates where the CCD falls 
ra_values_to_remove = [(343.55, 62.66, 63.36), (344.45, 62.68, 63.36), (344.44, 62.68, 63.36)]
tolerance = 1e-2  #RA precision tolerance

#Mask to exclude specified RA and DEC ranges
mask = np.ones(len(table), dtype=bool)

for ra, dec_min, dec_max in ra_values_to_remove:
    match_ra = np.isclose(table['RA'], ra, atol=tolerance)
    match_dec = (table['DEC'] >= dec_min) & (table['DEC'] <= dec_max)
    mask &= ~(match_ra & match_dec)

#Apply the mask
filtered_table = table[mask]
#Save filtered table
filtered_table.write(folder_path + "filtered_CCD_variables.vot", format="votable", overwrite=True)

print("Filtering complete! Saved as 'filtered_CCD_variables.vot'.")


