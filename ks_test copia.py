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

variable_stars_votable = "Simbad_variable_objects/variables_g_J_ws.vot"
non_variable_stars_votable = "Gaia_nonvariable_objects/nonvariable_g_J_ws.vot"
index_column = "ws_index"

# Read variability indices
variable_indices = read_variability_indices(variable_stars_votable, index_column)
non_variable_indices = read_variability_indices(non_variable_stars_votable, index_column)

# Perform the two-sample KS test
statistic, p_value = ks_2samp(variable_indices, non_variable_indices)

# Print results
print(f"KS Statistic: {statistic}")
print(f"P-Value: {p_value}")


# if p_value < 0.05:
#     print("The distributions of the variability index for variable and non-variable stars are significantly different.")
# else:
#     print("The distributions of the variability index for variable and non-variable stars are not significantly different.")


########## Figure ##########

# import numpy as np
# import matplotlib.pyplot as plt

# # Function to compute ECDF
# def compute_ecdf(data):
#     n = len(data)
#     sorted_data = np.sort(data)
#     ecdf = np.arange(1, n + 1) / n
#     return sorted_data, ecdf

# # Compute ECDF for variable and non-variable stars
# variable_sorted, variable_ecdf = compute_ecdf(variable_indices)
# non_variable_sorted, non_variable_ecdf = compute_ecdf(non_variable_indices)

# # Plot the ECDFs
# plt.figure(figsize=(8, 6))
# plt.plot(variable_sorted, variable_ecdf, label="Variable Stars", marker='.', linestyle='-', color='blue')
# plt.plot(non_variable_sorted, non_variable_ecdf, label="Non-Variable Stars", marker='.', linestyle='-', color='orange')

# # Customize the plot
# plt.title("Empirical Cumulative Distribution Functions (ECDF)")
# plt.xlabel("Variability Index")
# plt.ylabel("Cumulative Probability")
# plt.legend(loc="best")
# plt.grid(True)

# # Show the plot
# plt.show()
