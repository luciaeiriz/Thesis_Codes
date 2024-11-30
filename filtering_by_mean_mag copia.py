import numpy as np
import os
from astropy.table import Table
from astropy.io.votable import parse, from_table, writeto

# Paths
data_folder = 'Data_per_night'
output_folder = 'Filtered_Data'

# Create output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# List all files 
files = os.listdir(data_folder)

# Iterate over each file 
for target_file in files:
    if 'uJAVA' in target_file and target_file.endswith('.vot'):  # Check for 'gSDSS' in filename and .vot extension
        file_path = os.path.join(data_folder, target_file)

        votable = parse(file_path)
        table = votable.get_first_table().to_table()

        # Identify columns that contain 'mag' in header
        mag_columns = [col for col in table.colnames if 'mag_' in col.lower()]

        if mag_columns:
            # Extract magnitude data
            mag_matrix = np.array([table[col].data for col in mag_columns], dtype='float').T

            # Replace masked/empty values with NaN
            mag_matrix = np.where(np.ma.is_masked(mag_matrix), np.nan, mag_matrix)

            # Calculate mean magnitude for each row, considering only valid (non-NaN) values
            mean_magnitude = np.full(mag_matrix.shape[0], np.nan)
            for i in range(mag_matrix.shape[0]):
                valid_mags = mag_matrix[i][~np.isnan(mag_matrix[i])]
                if len(valid_mags) > 0:
                    mean_magnitude[i] = np.mean(valid_mags)

            # Filter rows where mean magnitude is x or above
            mask = mean_magnitude >= 14.4
            filtered_table = table[mask]

            # Create a new VOTable from the filtered data
            new_votable = from_table(filtered_table)

            # Write the filtered data to a new file in the output folder
            base_name = target_file.replace('_all_xcalibrated_clean_processed.vot', '')
            output_file = os.path.join(output_folder, f"{base_name}_filtered.vot")
            writeto(new_votable, output_file)

            print(f"Processed and filtered: {target_file}")
            print(f"Original rows: {len(table)}, Filtered rows: {len(filtered_table)}")
        else:
            print(f"No 'mag_' columns found in {target_file}")

print("Processing complete.")