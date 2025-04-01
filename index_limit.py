from astropy.table import Table
import os

filter = 'J0660' #change when necessary
index = 'sa04' #change when necessary

input_folder = 'SA04_Indices' #change when necessary
input_file = filter+'_'+index+'.vot'

output_file = 'filtered_'+ input_file
output_folder = 'Variable_objects/SA04_Variables' #change when necessary

input_path = os.path.join(input_folder, input_file)
output_path = os.path.join(output_folder, output_file)

# Load the VOTable
try:
    table = Table.read(input_path, format='votable')
except Exception as e:
    print(f"Error reading the VOTable: {e}")
    exit()

# Filter data
filtered_table = table[table.columns[3] >= 0.37] # above given limit #change when necessary
#filtered_table = table[(table.columns[3] < - 0.3) | (table.columns[3] > 0.35)] # below and above given limit #change when necessary

# Save the filtered table to a new VOTable
try:
    filtered_table.write(output_path, format='votable', overwrite=True)
    print(f"Filtered VOTable saved to {output_path}")
except Exception as e:
    print(f"Error writing the filtered VOTable: {e}")