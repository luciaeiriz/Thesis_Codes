from astropy.io.votable import parse, from_table, writeto
from astropy.table import Table
import numpy as np
import os

def calculate_welch_stetson(b, v, sigma_b, sigma_v):
    mask = ~np.isnan(b) & ~np.isnan(v) & ~np.isnan(sigma_b) & ~np.isnan(sigma_v) & (sigma_b > 0) & (sigma_v > 0)
    b, v, sigma_b, sigma_v = b[mask], v[mask], sigma_b[mask], sigma_v[mask]
    
    n = len(b)
    if n < 2:
        return np.nan, n, np.nan

    b_bar = np.sum(b / sigma_b**2) / np.sum(1 / sigma_b**2)
    v_bar = np.sum(v / sigma_v**2) / np.sum(1 / sigma_v**2)
    
    numerator = np.sum(((b - b_bar) / sigma_b) * ((v - v_bar) / sigma_v))
    denominator = np.sqrt(n * (n - 1))
    
    ws_index = numerator / denominator
    mean_mag = (b_bar + v_bar) / 2

    return ws_index, n, mean_mag

# Input and output folder setup
input_folder = 'Paired_data_by_filters'
output_folder = 'WS_Indices'
os.makedirs(output_folder, exist_ok=True)

# Process each VOTable in the input folder
for filename in os.listdir(input_folder):
    if filename.endswith('.vot'):
        input_file = os.path.join(input_folder, filename)
        try:
            votable = parse(input_file)
            data = votable.get_first_table().array
        except Exception as e:
            print(f"Error reading VOTable {filename}: {e}")
            continue

        # Initialize lists to store results
        ra_list, dec_list = [], []
        variability_indices, n_observations, mean_magnitudes = [], [], []

        # Iterate through each stellar object (row)
        for row in data:
            ra_list.append(row['RA_1'])
            dec_list.append(row['DEC_1'])
            
            # Extract magnitude and error data for both bands
            b_mags, v_mags, b_errors, v_errors = [], [], [], []
            for i in range(2, 67):
                if row[f'night_{i}_1'] == row[f'night_{i}_2']:  # Match observations by night
                    b_mags.append(row[f'mag_{i}_1'])
                    v_mags.append(row[f'mag_{i}_2'])
                    b_errors.append(row[f'mager_{i}_1'])
                    v_errors.append(row[f'mager_{i}_2'])
            
            # Convert to numpy arrays
            b_mags = np.ma.filled(np.ma.array(b_mags), fill_value=np.nan)
            v_mags = np.ma.filled(np.ma.array(v_mags), fill_value=np.nan)
            b_errors = np.ma.filled(np.ma.array(b_errors), fill_value=np.nan)
            v_errors = np.ma.filled(np.ma.array(v_errors), fill_value=np.nan)   
            
            # Calculate Welch-Stetson Variability Index
            index, n, mean_mag = calculate_welch_stetson(b_mags, v_mags, b_errors, v_errors)
            variability_indices.append(index)
            n_observations.append(n)
            mean_magnitudes.append(mean_mag)

        # Create a new table with the results
        new_table = Table(data=[ra_list, dec_list, n_observations, variability_indices, mean_magnitudes],
                          names=('RA', 'DEC', 'N', 'ws_index', 'mean_magnitude'))

        # Create a new VOTable from the table
        new_votable = from_table(new_table)

        # Generate the output file name
        base_name_without_ext = os.path.splitext(filename)[0]
        output_file_path = os.path.join(output_folder, f"{base_name_without_ext}_ws.vot")

        # Save the new VOTable
        try:
            writeto(new_votable, output_file_path)
            print(f"New VOTable saved as {output_file_path}")
        except Exception as e:
            print(f"Error writing VOTable {filename}: {e}")