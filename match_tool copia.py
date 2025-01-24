import os
from astropy.io.votable import parse, writeto
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import units as u
from astropy.table import Table

#Max error of 1.2 arcsec
max_error = 1.2 * u.arcsec

#Load VOTable. Change depnding on what is being processed!
gaia_votable = parse('Gaia_nonvariables.vot')
gaia_table = gaia_votable.get_first_table().to_table()

#Extract ra and dec from table being processed. Change name depending on table (Simbad/ Gaia)
gaia_coords = SkyCoord(ra=gaia_table['ra'], dec=gaia_table['dec'], unit=(u.degree, u.degree), frame='icrs')

#Output folder
output_folder = 'Gaia_nonvariable_objects'
os.makedirs(output_folder, exist_ok=True)
current_directory = os.getcwd()

# Loop through folder that end with Indices
for item in os.listdir(current_directory):
    if os.path.isdir(item) and item.endswith('_Indices'):
        print(f"Processing folder: {item}")
        
        for filename in os.listdir(item):
            if filename.endswith('.vot'):
                print(f"Processing file: {filename}")
                
                # Load file inside the folder
                chi_squared_votable = parse(os.path.join(item, filename))
                chi_squared_table = chi_squared_votable.get_first_table().to_table()

                # Extract coordinated for file being loaded 
                chi_squared_coords = SkyCoord(ra=chi_squared_table['RA'], dec=chi_squared_table['DEC'], unit=(u.degree, u.degree), frame='icrs')

                #Matches between both files
                idx_chi_squared, idx_gaia, d2d, _ = search_around_sky(chi_squared_coords, gaia_coords, max_error)

                #Creates new matched table
                matched_chi_squared = chi_squared_table[idx_chi_squared]
                matched_gaia = gaia_table[idx_gaia]

                #Rename columns in the Gaia table to avoid conflicts
                for colname in matched_gaia.colnames:
                    if colname in matched_chi_squared.colnames:
                        matched_gaia.rename_column(colname, f'gaia_{colname}')

                #Combine matched data into a new table
                result_table = Table(matched_chi_squared)
                for colname in matched_gaia.colnames:
                    result_table[colname] = matched_gaia[colname]

                # Add separation column
                result_table['separation'] = d2d.to(u.arcsec)

                # Save the matched data as a new VOTable
                output_filename = f"nonvariables_{filename}" #adds prefix to saved tables. Change depending on data being matched
                output_filepath = os.path.join(output_folder, output_filename)
                
                writeto(result_table, output_filepath)
                print(f"Saved matched table: {output_filepath}")
print("Matching complete!")

