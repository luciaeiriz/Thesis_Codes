import os
from astropy.io.votable import parse, writeto
from astropy.table import Table
from astropy.io.votable.tree import VOTableFile, Resource, Table as TableElement

# Path to the folder containing the original VOTables
folder_path = 'Chi_squared_indices' # Change when necessary

# Loop through each file in the folder
for file_name in os.listdir(folder_path):
    if file_name.endswith('.vot'):
        file_path = os.path.join(folder_path, file_name)
        
        votable = parse(file_path)
        table = votable.get_first_table()
        
        astropy_table = Table(table.array)
        
        if file_name.startswith('uJAVA'):
            threshold = 13
        else:  # For gSDSS, rSDSS, zSDSS, iSDSS, J0660
            threshold = 40
        
        filtered_table = astropy_table[astropy_table['N'] >= threshold] # Filter rows based on the threshold
        
        # Check if any rows remain after filtering
        if len(filtered_table) > 0:
            new_votable = VOTableFile() # Create a new filtered VOTable
            resource = Resource()
            new_votable.resources.append(resource)
            new_table = TableElement.from_table(new_votable, filtered_table)
            resource.tables.append(new_table)
            
            writeto(new_votable, file_path)
        else:
            print(f"No rows left after filtering for {file_name}. Skipping writing.")

print("Filtered VOTables successfully overwritten in the same folder.")
