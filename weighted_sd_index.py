import numpy as np
import os
from astropy.table import Table
from astropy.io.votable import parse, from_table, writeto

# Weighted Standard Deviation

# Paths
data_folder = 'Filtered_Data'
output_folder = 'Weighted_SD_Indices'

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

        # Weighted mean magnitude calculation
        weighted_mag = np.zeros(mag_matrix.shape[0])
        mean_magnitude = np.full(mag_matrix.shape[0], np.nan)
        for i in range(mag_matrix.shape[0]):
            # Extract valid magnitudes and errors (remove NaN)
            valid_mags = mag_matrix[i, ~np.isnan(mag_matrix[i]) & ~np.isnan(mager_matrix[i])]
            valid_mager = mager_matrix[i, ~np.isnan(mag_matrix[i]) & ~np.isnan(mager_matrix[i])]

            if len(valid_mags) > 0:  # Perform weighted mean calculation

                # Summation components for weighted mean
                numerator = np.sum(valid_mags / valid_mager ** 2)
                denominator = np.sum(1 / (valid_mager ** 2))

                # Weighted average
                weighted_mag[i] = numerator / denominator
                mean_magnitude[i] = np.nanmean(valid_mags)
            else:
                weighted_mag[i] = np.nan
                mean_magnitude[i] = np.nan

        # Weighted Standard Deviation based on the new formula
        weighted_sd = np.zeros(mag_matrix.shape[0])
        for i in range(mag_matrix.shape[0]):
            # Extract valid magnitudes and errors (remove NaN)
            valid_mask = ~np.isnan(mag_matrix[i]) & ~np.isnan(mager_matrix[i])
            valid_mags = mag_matrix[i, valid_mask]
            valid_mager = mager_matrix[i, valid_mask]

            # Calculate weight (w_i = 1 / sigma_i^2)
            w_i = 1 / valid_mager ** 2

            # Weighted mean calculation
            numerator = np.sum(valid_mags / valid_mager ** 2)
            denominator = np.sum(1 / (valid_mager ** 2))
            weighted_mag[i] = numerator / denominator

            if len(valid_mags) > 0:
                mean_magnitude[i] = np.nanmean(valid_mags)

                # Calculate components of the formula
                a = np.sum(w_i)  # Sum of weights
                b = np.sum(w_i) ** 2  # Square of the sum of weights
                c = np.sum(w_i ** 2)  # Sum of squared weights
                d = np.sum(w_i * (valid_mags - weighted_mag[i]) ** 2)  # Weighted sum of squared deviations

                # Calculate the weighted standard deviation using the formula
                weighted_sd[i] = np.sqrt((a * d) / (b - c))

            else:
                mean_magnitude[i] = np.nan
                weighted_sd[i] = np.nan

        # Saving data in new table
        ra_column = table['RA'].data
        dec_column = table['DEC'].data

        new_table = Table(data=[ra_column, dec_column, N, weighted_sd, mean_magnitude], names=('RA', 'DEC', 'N', 'weighted_standard_deviation', 'mean_magnitude'))
        new_votable = from_table(new_table)

        # Generate the output file name by removing '_filtered.vot' from the original name
        base_name = target_file.replace('.vot', '')
        output_file_path = os.path.join(output_folder, f"{base_name}_weightedsd.vot")

        # Save the new VOTable to the specified folder
        writeto(new_votable, output_file_path)

        print(f"New VOTable saved as {output_file_path}")

