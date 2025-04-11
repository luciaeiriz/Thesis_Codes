# Filter by number of nights
import os
import numpy as np
from astropy.io.votable import parse, writeto, from_table
from astropy.table import Table

# Paths
index_name = 'Stetson_K'
data_folder = index_name + '_Indices'
output_folder = data_folder

# Create output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# List all files
files = os.listdir(data_folder)

# Iterate over each file
for target_file in files:
    if target_file.endswith('.vot'):
        file_path = os.path.join(data_folder, target_file)

        # Parse the VOTable file
        votable = parse(file_path)
        table = votable.get_first_table().to_table()

        # Access the third column directly by index (number of nights)
        third_column = table.columns[2]  # Third column (index 2)

        # Check if the column data is numeric (for comparison)
        if np.issubdtype(third_column.dtype, np.number):
            if 'uJAVA' in target_file:
                threshold = 10  # Threshold of 10 nights for uJAVA
            elif 'J0660' in target_file or 'SDSS' in target_file:
                threshold = 25  # Threshold of 25 nights for all other fils
            else:
                threshold = None

            # Apply filtering if a threshold is set
            if threshold is not None:
                mask = third_column.data >= threshold
                filtered_table = table[mask]
            else:
                filtered_table = table  # No filtering if no condition met

            # Create a new VOTable from the filtered data
            new_votable = from_table(filtered_table)

            # Write the filtered data to a new file in the output folder
            base_name = target_file.replace('.vot', '')
            output_file = os.path.join(output_folder, f"{base_name}.vot")
            writeto(new_votable, output_file)

            print(f"Processed and filtered: {target_file}")
            print(f"Original rows: {len(table)}, Filtered rows: {len(filtered_table)}")
        else:
            print(f"Third column is not numeric in {target_file}")

print("Processing complete.")


# # filtering for paired data
# import os
# import numpy as np
# from astropy.io.votable import parse, writeto, from_table
# from astropy.table import Table

# index_name = 'Welch_Stetson'
# data_folder = index_name+'_Indices'
# output_folder = data_folder

# # Create output folder if it doesn't exist
# if not os.path.exists(output_folder):
#     os.makedirs(output_folder)

# # List all files 
# files = os.listdir(data_folder)

# # Iterate over each file 
# for target_file in files:
#     if '_' in target_file and target_file.endswith('.vot'):  # Check for 'uJAVA' in filename and .vot extension
#         file_path = os.path.join(data_folder, target_file)

#         # Parse the VOTable file
#         votable = parse(file_path)
#         table = votable.get_first_table().to_table()

#         # Access the third column directly by index
#         third_column = table.columns[2]  # Third column (index 2)

#         # Check if the column data is numeric (for comparison)
#         if np.issubdtype(third_column.dtype, np.number):
#             # Create a mask for rows where the third column value is >= 10
#             mask = third_column.data >= 25

#             # Filter the table using the mask
#             filtered_table = table[mask]

#             # Create a new VOTable from the filtered data
#             new_votable = from_table(filtered_table)

#             # Write the filtered data to a new file in the output folder
#             base_name = target_file.replace('.vot', '')
#             output_file = os.path.join(output_folder, f"{base_name}.vot")
#             writeto(new_votable, output_file)

#             print(f"Processed and filtered: {target_file}")
#             print(f"Original rows: {len(table)}, Filtered rows: {len(filtered_table)}")
#         else:
#             print(f"Third column is not numeric in {target_file}")