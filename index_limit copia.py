from astropy.table import Table
import os

input_folder = 'Stetson_J_Indices'
input_file = 'g_i_sJ.vot'
output_file = 'variable_g_i_sJ.vot'
output_folder = 'Variable_objects'

input_path = os.path.join(input_folder, input_file)
output_path = os.path.join(output_folder, output_file)

# Load the VOTable
try:
    table = Table.read(input_path, format='votable')
except Exception as e:
    print(f"Error reading the VOTable: {e}")
    exit()

# Filter data
#filtered_table = table[table.columns[3] >= 1.3]
filtered_table = table[(table.columns[3] < - 0.8) | (table.columns[3] > 0.8)]

# Save the filtered table to a new VOTable
try:
    filtered_table.write(output_path, format='votable', overwrite=True)
    print(f"Filtered VOTable saved to {output_path}")
except Exception as e:
    print(f"Error writing the filtered VOTable: {e}")