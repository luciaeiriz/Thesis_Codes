from astropy.table import Table
import numpy as np

# Define the path to the file
folder_path = 'Clean_Data'
file_name = 'CepOB3_rSDSS_NPh_aperphot_ststar_v4_part15_all_xcalibrated_clean.vot'
file_path = f"{folder_path}/{file_name}"

# Load the VOTable file using Astropy
data_table = Table.read(file_path, format='votable')

# Access row 5298 (zero-indexed as 5297)
row = data_table[4793]

# Extract all columns that contain 'Cal_mag_' in their names
cal_mag_columns = [col for col in data_table.colnames if 'Cal_mag_' in col]
magnitude = [row[col] if row[col] is not None else np.nan for col in cal_mag_columns]

# Extract all columns that contain 'Cal_mager_' in their names
cal_mager_columns = [col for col in data_table.colnames if 'Cal_mager_' in col]
magnitude_err = [row[col] if row[col] is not None else np.nan for col in cal_mager_columns]

# Convert the results to numpy arrays for easier handling
magnitudes = np.array(magnitude)
magnitude_errs = np.array(magnitude_err)

# Calculate minimum, maximum, and average, ignoring NaN values
mag_min = np.nanmin(magnitudes)
mag_max = np.nanmax(magnitudes)
mag_avg = np.nanmean(magnitudes)

mag_err_min = np.nanmin(magnitude_errs)
mag_err_max = np.nanmax(magnitude_errs)
mag_err_avg = np.nanmean(magnitude_errs)

# Display the results
print("Magnitude Min:", mag_min)
print("Magnitude Max:", mag_max)
print("Magnitude Avg:", mag_avg)
print()
print("Magnitude Error Min:", mag_err_min)
print("Magnitude Error Max:", mag_err_max)
print("Magnitude Error Avg:", mag_err_avg)



