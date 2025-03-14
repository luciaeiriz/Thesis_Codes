import os
import numpy as np
from astropy.io.votable import parse, from_table, writeto
from astropy.table import Table

data_folder = 'Stetson_Indices'
output_folder = 'Stetson_Indices'

# Create output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# List all files 
files = os.listdir(data_folder)

# Iterate over each file 
for target_file in files:
    if '_' in target_file and target_file.endswith('.vot'):  # Check for 'filter' in filename and .vot extension
        file_path = os.path.join(data_folder, target_file) 

        votable = parse(file_path)
        table = votable.get_first_table().to_table()

        # Access the third column (index column)
        third_column = table.columns[2]

        if np.issubdtype(third_column.dtype, np.number):
            
            mask = third_column.data >= 25 # 10 for uJAVA, 25 for all other files
            filtered_table = table[mask] # Filter the table using the mask

            new_votable = from_table(filtered_table) # Create a new VOTable from the filtered data

            # Write the filtered data to a new file in the output folder
            base_name = target_file.replace('.vot', '')
            output_file = os.path.join(output_folder, f"{base_name}_filtered.vot")
            writeto(new_votable, output_file)

            print(f"Processed and filtered: {target_file}")
            print(f"Original rows: {len(table)}, Filtered rows: {len(filtered_table)}")
        else:
            print(f"Third column is not numeric in {target_file}")

print("Processing complete.")