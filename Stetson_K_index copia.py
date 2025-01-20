import numpy as np
import os
from astropy.table import Table
from astropy.io.votable import parse, from_table, writeto

# Stetson K Index 

# Paths
data_folder = 'Filtered_Data'
output_folder = 'Stetson_K_Indices'

# List all files 
files = os.listdir(data_folder)

# Iterate over each file 
for target_file in files:
    if target_file.endswith('.vot'):  # Check for specific file 
        file_path = os.path.join(data_folder, target_file)

        votable = parse(file_path)
        table = votable.get_first_table().to_table()

        # Magnitude data
        magnitude = []
        magnitude_err = []

        # Identify columns that contain '_mag_' and '_mager_' in header
        for col in table.colnames:
            if 'mag_' in col:
                magnitude.append(table[col].data)
            elif 'mager_' in col:
                magnitude_err.append(table[col].data)

        # Convert lists to numpy arrays, replacing empty values with NaN
        mag_matrix = np.array(magnitude, dtype='float').T  # Transpose for objects in rows
        mager_matrix = np.array(magnitude_err, dtype='float').T  # Transpose for objects in rows

        # Replace masked/empty values with NaN
        mag_matrix = np.where(np.ma.is_masked(mag_matrix), np.nan, mag_matrix)
        mager_matrix = np.where(np.ma.is_masked(mager_matrix), np.nan, mager_matrix)

        # Number of nights for each star, ignoring those that are empty
        N = np.sum(~np.isnan(mag_matrix), axis=1)

        # Stetson K Index calculation
        sk_index = np.full(mag_matrix.shape[0], np.nan)  # Initialize Stetson K values
        mean_magnitude = np.full(mag_matrix.shape[0], np.nan)  # Initialize mean magnitudes
        
        for i in range(mag_matrix.shape[0]):
            # Extract valid magnitudes and errors (remove NaN)
            valid_mask = ~np.isnan(mag_matrix[i]) & ~np.isnan(mager_matrix[i])
            valid_mags = mag_matrix[i, valid_mask]
            valid_mager = mager_matrix[i, valid_mask]

            if len(valid_mags) > 1: 
                
                weights = 1 / valid_mager**2  # Weights
                weighted_mean = np.sum(valid_mags * weights) / np.sum(weights)

                delta_i = np.sqrt(len(valid_mags) / (len(valid_mags) - 1)) * (valid_mags - weighted_mean) / valid_mager

                numerator = np.sum(np.abs(delta_i)) / len(valid_mags)
                denominator = np.sqrt(np.sum(delta_i**2) / len(valid_mags))
                sk_index[i] = numerator / denominator

                # Store mean magnitude
                mean_magnitude[i] = weighted_mean
            else:
                # Not enough data points
                sk_index[i] = np.nan
                mean_magnitude[i] = np.nan

        # Saving data in new table
        ra_column = table['RA'].data
        dec_column = table['DEC'].data

        new_table = Table(data=[ra_column, dec_column, N, sk_index, mean_magnitude], names=('RA', 'DEC', 'N', 'Stetson_K', 'mean_magnitude'))
        new_votable = from_table(new_table)

        # Generate the output file name by removing '_filtered.vot' from the original name
        base_name = target_file.replace('.vot', '')
        output_file_path = os.path.join(output_folder, f"{base_name}_sk.vot")

        # Save the new VOTable to the specified folder
        writeto(new_votable, output_file_path)

        print(f"New VOTable saved as {output_file_path}")
