import numpy as np
import os
from astropy.table import Table
from astropy.io.votable import parse, from_table, writeto

# Lag 1 index

# Paths
data_folder = 'Filtered_Data'
output_folder = 'Lag1_Indices'

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

        # Initialize arrays to store results
        l1_index = np.zeros(mag_matrix.shape[0])
        mean_magnitude = np.full(mag_matrix.shape[0], np.nan)
        weighted_mag = np.full(mag_matrix.shape[0], np.nan)

        for i in range(mag_matrix.shape[0]):
            # Extract valid magnitudes and errors (remove NaN)
            valid_mags = mag_matrix[i, ~np.isnan(mag_matrix[i]) & ~np.isnan(mager_matrix[i])]
            valid_mager = mager_matrix[i, ~np.isnan(mag_matrix[i]) & ~np.isnan(mager_matrix[i])]

            if len(valid_mags) > 1:  # Perform calculation only if there are at least 2 values
                # Calculate mean magnitude
                mean_magnitude[i] = np.nanmean(valid_mags)

                # Calculate weighted mean magnitude
                numerator_weighted_mean = np.sum(valid_mags / valid_mager**2)
                denominator_weighted_mean = np.sum(1 / valid_mager**2)
                weighted_mag[i] = numerator_weighted_mean / denominator_weighted_mean

                # Lag-1 index calculation
                # Numerator: sum of (m_i - weighted_mean) * (m_{i+1} - weighted_mean)
                numerator = np.sum((valid_mags[:-1] - weighted_mag[i]) * (valid_mags[1:] - weighted_mag[i]))
                # Denominator: sum of (m_i - weighted_mean)^2
                denominator = np.sum((valid_mags - weighted_mag[i]) ** 2)

                # Calculate l1 index, checking if the denominator is non-zero
                l1_index[i] = numerator / denominator if denominator != 0 else np.nan
            else:
                # If there aren't enough valid data points, store NaN
                mean_magnitude[i] = np.nan
                weighted_mag[i] = np.nan
                l1_index[i] = np.nan

        # Saving data in new table
        ra_column = table['RA'].data
        dec_column = table['DEC'].data

        new_table = Table(data=[ra_column, dec_column, N, l1_index, mean_magnitude],
                          names=('RA', 'DEC', 'N', 'l1_index', 'mean_magnitude'))
        new_votable = from_table(new_table)

        # Generate the output file name by removing '_filtered.vot' from the original name
        base_name = target_file.replace('.vot', '')
        output_file_path = os.path.join(output_folder, f"{base_name}_l1.vot")

        # Save the new VOTable to the specified folder
        writeto(new_votable, output_file_path)
        print(f"New VOTable saved as {output_file_path}")
