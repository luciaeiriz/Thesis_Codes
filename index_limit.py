from astropy.table import Table
import os

filter = 'r_z' # Change when necessary
index = 'stetson' # Change when necessary

input_folder = 'Stetson_Indices' # Change when necessary
input_file = filter+'_'+index+'.vot'

output_file = 'filtered_'+ input_file
output_folder = 'Variable_objects'

input_path = os.path.join(input_folder, input_file)
output_path = os.path.join(output_folder, output_file)

try:
    table = Table.read(input_path, format='votable')
except Exception as e:
    print(f"Error reading the VOTable: {e}")
    exit()

# Filter data above given limit
filtered_table = table[table.columns[3] >= 2] #change when necessary
# Filter data below and above given limit
#filtered_table = table[(table.columns[3] < - 0.5) | (table.columns[3] > 0.8)] #change when necessary

# Save the filtered table to a new VOTable
try:
    filtered_table.write(output_path, format='votable', overwrite=True)
    print(f"Filtered VOTable saved to {output_path}")
except Exception as e:
    print(f"Error writing the filtered VOTable: {e}")