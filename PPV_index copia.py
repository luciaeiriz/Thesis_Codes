import numpy as np
import os
from astropy.table import Table
from astropy.io.votable import parse, from_table, writeto

# Peak to peak variability

# Paths
data_folder = 'Filtered_Data'
output_folder = 'PPV_Indices'

# List all files 
files = os.listdir(data_folder)

# Iterate over each file 
for target_file in files:
    if target_file.endswith('_filtered.vot'):  # Check for specific file 
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

        # MAD
        mean_magnitude = np.full(mag_matrix.shape[0], np.nan)  # Array to store mean magnitudes
        peak_to_peak = np.zeros(mag_matrix.shape[0])
        for i in range(mag_matrix.shape[0]):
            # Extract valid magnitudes and errors (remove NaN)
            valid_mags = mag_matrix[i, ~np.isnan(mag_matrix[i]) & ~np.isnan(mager_matrix[i])]
            valid_mager = mager_matrix[i, ~np.isnan(mag_matrix[i]) & ~np.isnan(mager_matrix[i])]

            if len(valid_mags) > 0:  # Perform calculation only if there are valid values
                # Calculate mean magnitude
                mean_magnitude[i] = np.nanmean(valid_mags)

                max_bound = np.max(valid_mags - valid_mager)
                min_bound = np.min(valid_mags + valid_mager)

                numerator = max_bound - min_bound
                denominator = max_bound + min_bound
                peak_to_peak[i] = numerator / denominator if denominator != 0 else np.nan
            else:
                # If no valid data, store NaN
                mean_magnitude[i] = np.nan
                peak_to_peak[i] = np.nan

        # Saving data in new table
        ra_column = table['RA'].data
        dec_column = table['DEC'].data

        new_table = Table(data=[ra_column, dec_column, N, peak_to_peak, mean_magnitude], names=('RA', 'DEC', 'N', 'peak_to_peak_variability', 'mean_magnitude'))
        new_votable = from_table(new_table)

        # Generate the output file name by removing '_all_xcalibrated_clean.vot' from the original name
        base_name = target_file.replace('_filtered.vot', '')
        output_file_path = os.path.join(output_folder, f"{base_name}_ppv.vot")

        # Save the new VOTable to the specified folder
        writeto(new_votable, output_file_path)

        print(f"New VOTable saved as {output_file_path}")