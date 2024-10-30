import numpy as np
import os
from astropy.table import Table
from astropy.io.votable import parse, from_table, writeto

# Chi Squared Index

# Paths
data_folder = 'Clean_Data'
output_folder = 'Chi_squared_Indices'

# List all files 
files = os.listdir(data_folder)

# Iterate over each file 
for target_file in files:
    if target_file.endswith('_all_xcalibrated_clean.vot'):  # Check for specific file 
        file_path = os.path.join(data_folder, target_file)

        votable = parse(file_path)
        table = votable.get_first_table().to_table()

        # Magnitude data
        magnitude = []
        magnitude_err = []

        # Identify columns that contain '_mag_' and '_mager_' in header
        for col in table.colnames:
            if '_mag_' in col:
                magnitude.append(table[col].data)
            elif '_mager_' in col:
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
        for i in range(mag_matrix.shape[0]):
            # Extract valid magnitudes and errors (remove NaN)
            valid_mags = mag_matrix[i, ~np.isnan(mag_matrix[i]) & ~np.isnan(mager_matrix[i])]
            valid_mager = mager_matrix[i, ~np.isnan(mag_matrix[i]) & ~np.isnan(mager_matrix[i])]

            if len(valid_mags) > 0:  # Perform calculation only if there are valid values
                # Summation components
                numerator = np.sum(valid_mags / valid_mager ** 2)
                denominator = np.sum(1 / (valid_mager ** 2))

                # Weighted average
                weighted_mag[i] = numerator / denominator
            else:
                # If no valid data, store NaN
                weighted_mag[i] = np.nan

        # Chi-squared
        chi_squared = np.zeros(mag_matrix.shape[0])
        for i in range(mag_matrix.shape[0]):
            # Extract valid magnitudes and errors (remove NaN)
            valid_mags = mag_matrix[i, ~np.isnan(mag_matrix[i]) & ~np.isnan(mager_matrix[i])]
            valid_mager = mager_matrix[i, ~np.isnan(mag_matrix[i]) & ~np.isnan(mager_matrix[i])]

            if len(valid_mags) > 0:  # Perform calculation only if there are valid values
                # Calculate chi-squared
                chi_squared[i] = np.sum(((valid_mags - weighted_mag[i]) ** 2) / (valid_mager ** 2))
                # Divide chi-squared by the number of valid points (N)
                chi_squared[i] /= N[i]  # Only divide if N > 0
            else:
                # If no valid data, store NaN
                chi_squared[i] = np.nan

        # Saving data in new table
        ra_column = table['RA'].data
        dec_column = table['DEC'].data

        new_table = Table(data=[ra_column, dec_column, N, chi_squared], names=('RA', 'DEC', 'N', 'chi_squared'))
        new_votable = from_table(new_table)

        # Generate the output file name by removing '_all_xcalibrated_clean.vot' from the original name
        base_name = target_file.replace('_all_xcalibrated_clean.vot', '')
        output_file_path = os.path.join(output_folder, f"{base_name}_chisquared.vot")

        # Save the new VOTable to the specified folder
        writeto(new_votable, output_file_path)

        print(f"New VOTable saved as {output_file_path}")
