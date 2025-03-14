import os
import pandas as pd
import numpy as np
from astropy.io.votable import parse_single_table, from_table, writeto
from astropy.table import Table

def process_vot_file(file_path):
    table = parse_single_table(file_path)
    
    if not isinstance(table, Table):
        table = table.to_table()
    
    df = table.to_pandas()

    result_df = df.iloc[:, :6].copy() # Ignore first 6 columns
    remaining_df = df.iloc[:, 6:] # Process the remaining columns

    # Identify columns with 'JD', 'mag', and 'mager' in their headers
    jd_columns = [col for col in remaining_df.columns if 'JD' in col]
    mag_columns = [col for col in remaining_df.columns if 'mag' in col.lower() and 'mager' not in col.lower()]
    mager_columns = [col for col in remaining_df.columns if 'mager' in col.lower()]

    def process_star(row):
        nights = {}
        for jd_col, mag_col, mager_col in zip(jd_columns, mag_columns, mager_columns):
            jd = row[jd_col]
            mag = row[mag_col]
            mager = row[mager_col]
            
            if pd.notna(jd) and pd.notna(mag) and pd.notna(mager):
                night = int(np.floor(jd))
                if night not in nights:
                    nights[night] = {'mags': [], 'magers': []}
                nights[night]['mags'].append(mag)
                nights[night]['magers'].append(mager)
        
        # Calculate weighted mean magnitude and weighted mean error for each night
        return {
            night: {
                'avg_mag': np.average(data['mags'], weights=1/np.array(data['magers'])**2),
                'avg_mager': np.sqrt(1 / np.sum(1/np.array(data['magers'])**2))
            }
            for night, data in nights.items()
        }

    processed_data = remaining_df.apply(process_star, axis=1)

    all_nights = sorted(set.union(*[set(star.keys()) for star in processed_data]))

    new_data = {} # Prepare data for new columns
    for i, night in enumerate(all_nights, start=2):
        new_data[f'night_{i}'] = [night] * len(processed_data)
        new_data[f'mag_{i}'] = processed_data.apply(lambda x: x.get(night, {}).get('avg_mag', np.nan))
        new_data[f'mager_{i}'] = processed_data.apply(lambda x: x.get(night, {}).get('avg_mager', np.nan))

    new_df = pd.DataFrame(new_data)
    result_df = pd.concat([result_df, new_df], axis=1)

    return result_df

# Folder path
folder_path = 'Clean_Data'
output_path = 'Data_per_night'

if not os.path.exists(output_path): # Ensure folder path exists
    os.makedirs(output_path)

# Loop over all .vot files in the 'Clean_Data' folder
for file_name in os.listdir(folder_path):
    if file_name.endswith('.vot'):
        file_path = os.path.join(folder_path, file_name)
        
        result = process_vot_file(file_path)
        astropy_table = Table.from_pandas(result)

        votable = from_table(astropy_table) # Create new VOTable

        output_file_name = file_name.replace('.vot', '.vot')
        output_file_path = os.path.join(output_path, output_file_name)

        # Save new file
        writeto(votable, output_file_path)
        print(f"Processing complete for {file_name}. Results saved to '{output_file_path}'")


