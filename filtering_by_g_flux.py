from astropy.table import Table


table = Table.read('simbad_variables.vot', format='votable')

# Filter out rows where gmeanmag > 19.5
filtered_table = table[table['FLUX_G'] <= 19.5]

# Save to a new VOTable file
filtered_table.write('filtered_simbad_variables.vot', format='votable', overwrite=True)