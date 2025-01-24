# #KS tets performed on all files at once
# import os
# import numpy as np
# from astropy.io.votable import parse
# from scipy.stats import ks_2samp
# import pandas as pd

# def read_variability_indices(votable_path):
#     votable = parse(votable_path)
#     table = votable.get_first_table()
#     data = table.array
#     return data[data.dtype.names[3]]

# def perform_ks_test(variable_file, nonvariable_file):
#     variable_indices = read_variability_indices(variable_file)
#     non_variable_indices = read_variability_indices(nonvariable_file)
#     statistic, p_value = ks_2samp(variable_indices, non_variable_indices)
#     return statistic, p_value

# variable_folder = 'Simbad_variable_objects'
# nonvariable_folder = 'Gaia_nonvariable_objects'

# results = []


# for filename in os.listdir(variable_folder):
#     if filename.startswith('variables_') and filename.endswith('.vot'):
#         identifier = filename.split('variables_')[1].split('.vot')[0]
        
#         variable_file = os.path.join(variable_folder, filename)
#         nonvariable_file = os.path.join(nonvariable_folder, f'nonvariables_{identifier}.vot')
        
#         if os.path.exists(nonvariable_file):
#             try:
#                 statistic, p_value = perform_ks_test(variable_file, nonvariable_file)
#                 results.append({
#                     'identifier': identifier,
#                     'ks_statistic': statistic,
#                     'p_value': p_value
#                 })
#             except Exception as e:
#                 print(f"Error processing {identifier}: {e}")

# results_df = pd.DataFrame(results)

# results_df.to_csv('ks_test_results.csv', index=False)


#KS test performed for two given files
import numpy as np
from astropy.io.votable import parse
from scipy.stats import ks_2samp

# Read variability indices from VOTable
def read_variability_indices(votable_path, index_column_name):
    votable = parse(votable_path)
    table = votable.get_first_table()
    data = table.array
    # Extract the column with variability index values
    return data[index_column_name].data 

variable_stars_votable = "Simbad_variable_objects/variables_J0660_chisquared.vot"
non_variable_stars_votable = "Gaia_nonvariable_objects/nonvariables_J0660_chisquared.vot"
index_column = "chi_squared"

# Read variability indices
variable_indices = read_variability_indices(variable_stars_votable, index_column)
non_variable_indices = read_variability_indices(non_variable_stars_votable, index_column)

# Perform the two-sample KS test
statistic, p_value = ks_2samp(variable_indices, non_variable_indices)

# Print results
print(f"KS Statistic: {statistic}")
print(f"P-Value: {p_value}")


if p_value < 0.05:
    print("The distributions of the variability index for variable and non-variable stars are significantly different.")
else:
    print("The distributions of the variability index for variable and non-variable stars are not significantly different.")


######### Figure ##########

# import numpy as np
# import matplotlib.pyplot as plt

# # Function for ECDF
# def compute_ecdf(data):
#     n = len(data)
#     sorted_data = np.sort(data)
#     ecdf = np.arange(1, n + 1) / n
#     return sorted_data, ecdf


# variable_sorted, variable_ecdf = compute_ecdf(variable_indices)
# non_variable_sorted, non_variable_ecdf = compute_ecdf(non_variable_indices)


# plt.figure(figsize=(8, 6))
# plt.plot(variable_sorted, variable_ecdf, label="Variable Stars", marker='.', linestyle='-', color='blue')
# plt.plot(non_variable_sorted, non_variable_ecdf, label="Non-Variable Stars", marker='.', linestyle='-', color='orange')

# plt.xscale('log')
# plt.title("Empirical Cumulative Distribution Functions (ECDF)")
# plt.xlabel("Chi Squared Index")
# plt.ylabel("Cumulative Probability")
# plt.legend(loc="best")
# plt.grid(True)

# plt.show()
