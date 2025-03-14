#KS tets performed on all files at once and produces a csv file with output 
import os
import numpy as np
from astropy.io.votable import parse
from scipy.stats import ks_2samp
import pandas as pd

def read_variability_indices(votable_path):
    votable = parse(votable_path)
    table = votable.get_first_table()
    data = table.array
    return data[data.dtype.names[3]]

def perform_ks_test(variable_file, nonvariable_file):
    variable_indices = read_variability_indices(variable_file)
    non_variable_indices = read_variability_indices(nonvariable_file)
    statistic, p_value = ks_2samp(variable_indices, non_variable_indices)
    return statistic, p_value

variable_folder = 'Simbad_variable_objects'
nonvariable_folder = 'Gaia_nonvariable_objects'

results = []


for filename in os.listdir(variable_folder):
    if filename.startswith('variables_') and filename.endswith('.vot'):
        identifier = filename.split('variables_')[1].split('.vot')[0]
        
        variable_file = os.path.join(variable_folder, filename)
        nonvariable_file = os.path.join(nonvariable_folder, f'nonvariables_{identifier}.vot')
        
        if os.path.exists(nonvariable_file):
            try:
                statistic, p_value = perform_ks_test(variable_file, nonvariable_file)
                results.append({
                    'identifier': identifier,
                    'ks_statistic': statistic,
                    'p_value': p_value
                })
            except Exception as e:
                print(f"Error processing {identifier}: {e}")

results_df = pd.DataFrame(results)

results_df.to_csv('ks_test_results.csv', index=False)