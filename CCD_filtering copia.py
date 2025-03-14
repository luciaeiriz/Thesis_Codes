from astropy.table import Table
import numpy as np 

table = Table.read("Variable_objects/Chi_Squared_Variables/variables_chisquared.vot", format="votable")

RA_CCD = [343.36, 343.55] # RA coordinates of the CCD 
DEC_min, DEC_max = 62.27, 63.38 # DEC coordinates range of the CCD
tolerance = 0.01  # maximum error deviation of 1.2 arc seconds

mask = ~((np.any([np.isclose(table['RA'], RA, atol=tolerance) for RA in RA_CCD], axis=0)) &
         (table['DEC'] >= DEC_min) & (table['DEC'] <= DEC_max))
filtered_table = table[mask]

filtered_table.write("filtered_file.vot", format="votable", overwrite=True) # save file 


