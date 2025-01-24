import os
from astropy.io.votable import parse, writeto
from astropy.table import Table
from astropy.io.votable.tree import VOTableFile, Resource, Table as TableElement

# Path to the folder containing the original VOTables
folder_path = 'Chi_squared_indices'

# Loop through each file in the folder
for file_name in os.listdir(folder_path):
    if file_name.endswith('.vot'):
        file_path = os.path.join(folder_path, file_name)
        
        # Parse the VOTable
        votable = parse(file_path)
        table = votable.get_first_table()
        
        # Convert to an Astropy Table for easier row filtering
        astropy_table = Table(table.array)
        
        # Determine the threshold based on file name prefix
        if file_name.startswith('uJAVA'):
            threshold = 13
        else:  # For gSDSS, rSDSS, zSDSS, iSDSS, J0660
            threshold = 40
        
        # Filter rows based on the threshold
        filtered_table = astropy_table[astropy_table['N'] >= threshold]
        
        # Check if any rows remain after filtering
        if len(filtered_table) > 0:
            # Create a new VOTable
            new_votable = VOTableFile()
            resource = Resource()
            new_votable.resources.append(resource)
            # Create a TableElement from the filtered table
            new_table = TableElement.from_table(new_votable, filtered_table)
            resource.tables.append(new_table)
            
            # Overwrite the original VOTable
            writeto(new_votable, file_path)
        else:
            print(f"No rows left after filtering for {file_name}. Skipping writing.")

print("Filtered VOTables successfully overwritten in the same folder.")
