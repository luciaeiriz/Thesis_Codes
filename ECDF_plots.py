import numpy as np
import matplotlib.pyplot as plt
from astropy.io.votable import parse
from scipy.stats import ks_2samp

# Function to read variability indices from VOTable
def read_variability_indices(votable_path, index_column_name):
    votable = parse(votable_path)
    table = votable.get_first_table()
    data = table.array[index_column_name].data
    
    # Filter out NaN and zero values
    filtered_data = data[~np.isnan(data) & (data > 0)]
    return filtered_data

# Function to compute ECDF
def compute_ecdf(data):
    n = len(data)
    sorted_data = np.sort(data)
    ecdf = np.arange(1, n + 1) / n
    return sorted_data, ecdf

# # 6 files
# x = "sa04" # Change when necessary
# file_pairs = [
#     ("Simbad_variable_objects/variables_gSDSS_"+x+".vot", "Gaia_nonvariable_objects/nonvariables_gSDSS_"+x+".vot"),
#     ("Simbad_variable_objects/variables_iSDSS_"+x+".vot", "Gaia_nonvariable_objects/nonvariables_iSDSS_"+x+".vot"),
#     ("Simbad_variable_objects/variables_zSDSS_"+x+".vot", "Gaia_nonvariable_objects/nonvariables_zSDSS_"+x+".vot"),
#     ("Simbad_variable_objects/variables_rSDSS_"+x+".vot", "Gaia_nonvariable_objects/nonvariables_rSDSS_"+x+".vot"),
#     ("Simbad_variable_objects/variables_J0660_"+x+".vot", "Gaia_nonvariable_objects/nonvariables_J0660_"+x+".vot"),
#     ("Simbad_variable_objects/variables_uJAVA_"+x+".vot", "Gaia_nonvariable_objects/nonvariables_uJAVA_"+x+".vot")
# ]

# titles = ["gSDSS", "iSDSS", "zSDSS", "rSDSS", "J0660", "uJAVA"]

# index_column = "Varia_night"  # Change when necessary

# # Create subplots
# fig, axes = plt.subplots(2, 3, figsize=(14, 8))
# axes = axes.flatten()

# # Iterate over file pairs and plot ECDFs
# for i, (var_file, nonvar_file) in enumerate(file_pairs):
#     variable_indices = read_variability_indices(var_file, index_column)
#     non_variable_indices = read_variability_indices(nonvar_file, index_column)
    
#     # Compute ECDFs
#     variable_sorted, variable_ecdf = compute_ecdf(variable_indices)
#     non_variable_sorted, non_variable_ecdf = compute_ecdf(non_variable_indices)
    
#     # Perform KS Test
#     statistic, p_value = ks_2samp(variable_indices, non_variable_indices)
    
#     # Extract title from filename
#     title = titles[i]
    
#     # Plot ECDFs
#     axes[i].plot(variable_sorted, variable_ecdf, label="Variable Stars", marker='.', linestyle='-', color='blue')
#     axes[i].plot(non_variable_sorted, non_variable_ecdf, label="Non-Variable Stars", marker='.', linestyle='-', color='red')
#     axes[i].legend(loc="upper left", fontsize=10)
#     axes[i].set_xscale('log')
#     axes[i].set_title(f"{title}\nKS p-value: {p_value:.3e}", fontsize = 11)
#     axes[i].set_xlabel("SA 04", fontsize = 13) # Change when necessary
#     axes[i].set_ylabel("Cumulative Probability", fontsize = 13)
#     axes[i].grid(True)
# plt.tight_layout()
# plt.savefig('Figures/ECDF_'+x)


# 10 files
x = "ws" # Change when necessary
file_pairs = [
    ("Simbad_variable_objects/variables_J_z_"+x+".vot", "Gaia_nonvariable_objects/nonvariables_J_z_"+x+".vot"),
    ("Simbad_variable_objects/variables_i_z_"+x+".vot", "Gaia_nonvariable_objects/nonvariables_i_z_"+x+".vot"),
    ("Simbad_variable_objects/variables_r_z_"+x+".vot", "Gaia_nonvariable_objects/nonvariables_r_z_"+x+".vot"),
    ("Simbad_variable_objects/variables_i_J_"+x+".vot", "Gaia_nonvariable_objects/nonvariables_i_J_"+x+".vot"),
    ("Simbad_variable_objects/variables_g_r_"+x+".vot", "Gaia_nonvariable_objects/nonvariables_g_r_"+x+".vot"),
    ("Simbad_variable_objects/variables_i_r_"+x+".vot", "Gaia_nonvariable_objects/nonvariables_i_r_"+x+".vot"),
    ("Simbad_variable_objects/variables_g_J_"+x+".vot", "Gaia_nonvariable_objects/nonvariables_g_J_"+x+".vot"),
    ("Simbad_variable_objects/variables_g_i_"+x+".vot", "Gaia_nonvariable_objects/nonvariables_g_i_"+x+".vot"),
    ("Simbad_variable_objects/variables_g_z_"+x+".vot", "Gaia_nonvariable_objects/nonvariables_g_z_"+x+".vot"),
    ("Simbad_variable_objects/variables_J_r_"+x+".vot", "Gaia_nonvariable_objects/nonvariables_J_r_"+x+".vot")
]

titles = ["J0660 & zSDSS", "iSDSS & zSDSS", "rSDSS & zSDSS", "iSDSS & J0660", "gSDSS & rSDSS", "iSDSS & rSDSS",
        "gSDSS & J0660", "gSDSS & iSDSS", "gSDSS & zSDSS", "J0660 & rSDSS"]

index_column = "Welch_Stetson"  # Change when necessary

# Create subplots
fig, axes = plt.subplots(2, 5, figsize=(19, 8))
axes = axes.flatten()

# Iterate over file pairs and plot ECDFs
for i, (var_file, nonvar_file) in enumerate(file_pairs):
    variable_indices = read_variability_indices(var_file, index_column)
    non_variable_indices = read_variability_indices(nonvar_file, index_column)
    
    # Compute ECDFs
    variable_sorted, variable_ecdf = compute_ecdf(variable_indices)
    non_variable_sorted, non_variable_ecdf = compute_ecdf(non_variable_indices)
    
    # Perform KS Test
    statistic, p_value = ks_2samp(variable_indices, non_variable_indices)
    
    # Extract title from filename
    title = titles[i]
    
    # Plot ECDFs
    axes[i].plot(variable_sorted, variable_ecdf, label="Variable Stars", marker='.', linestyle='-', color='blue')
    axes[i].plot(non_variable_sorted, non_variable_ecdf, label="Non-Variable Stars", marker='.', linestyle='-', color='red')
    axes[i].legend(loc="upper left", fontsize=8)
    axes[i].set_xscale('log')
    axes[i].set_title(f"{title}\nKS p-value: {p_value:.3e}", fontsize = 11)
    axes[i].set_xlabel("Welch Stetson", fontsize = 13) # Change when necessary
    axes[i].set_ylabel("Cumulative Probability", fontsize = 13)
    axes[i].grid(True)
plt.tight_layout()
plt.savefig('Figures/ECDF_'+x)


## ECDF subplots (?)
# import numpy as np
# import matplotlib.pyplot as plt
# from astropy.io.votable import parse
# from scipy.stats import ks_2samp

# # Function to read variability indices from VOTable
# def read_variability_indices(votable_path, index_column_name):
#     votable = parse(votable_path)
#     table = votable.get_first_table()
#     data = table.array
#     return data[index_column_name].data 

# # Function to compute ECDF
# def compute_ecdf(data):
#     n = len(data)
#     sorted_data = np.sort(data)
#     ecdf = np.arange(1, n + 1) / n
#     return sorted_data, ecdf

# # Define file pairs (6 categories, each with variable & non-variable file)
# file_pairs = [
#     ("Simbad_variable_objects/variables_rSDSS_chisquared.vot", "Gaia_nonvariable_objects/nonvariables_rSDSS_chisquared.vot"),
#     ("Simbad_variable_objects/variables_iSDSS_chisquared.vot", "Gaia_nonvariable_objects/nonvariables_iSDSS_chisquared.vot"),
#     ("Simbad_variable_objects/variables_zSDSS_chisquared.vot", "Gaia_nonvariable_objects/nonvariables_zSDSS_chisquared.vot"),
#     ("Simbad_variable_objects/variables_gSDSS_chisquared.vot", "Gaia_nonvariable_objects/nonvariables_gSDSS_chisquared.vot"),
#     ("Simbad_variable_objects/variables_J0660_chisquared.vot", "Gaia_nonvariable_objects/nonvariables_J0660_chisquared.vot"),
#     ("Simbad_variable_objects/variables_uJAVA_chisquared.vot", "Gaia_nonvariable_objects/nonvariables_uJAVA_chisquared.vot")
# ]

# index_column = "chi_squared"

# # Create subplots
# fig, axes = plt.subplots(2, 3, figsize=(12, 8))
# axes = axes.flatten()

# # Iterate over file pairs and plot ECDFs
# for i, (var_file, nonvar_file) in enumerate(file_pairs):
#     variable_indices = read_variability_indices(var_file, index_column)
#     non_variable_indices = read_variability_indices(nonvar_file, index_column)
    
#     # Compute ECDFs
#     variable_sorted, variable_ecdf = compute_ecdf(variable_indices)
#     non_variable_sorted, non_variable_ecdf = compute_ecdf(non_variable_indices)
    
#     # Perform KS Test
#     statistic, p_value = ks_2samp(variable_indices, non_variable_indices)
    
#     # Extract title from filename
#     title = var_file.split("_")[3]
    
#     # Plot ECDFs
#     axes[i].plot(variable_sorted, variable_ecdf, label="Variable Stars", marker='.', linestyle='-', color='blue')
#     axes[i].plot(non_variable_sorted, non_variable_ecdf, label="Non-Variable Stars", marker='.', linestyle='-', color='orange')
    
#     axes[i].set_xscale('log')
#     axes[i].set_title(f"{title}\nKS p-value: {p_value:.3e}")
#     axes[i].set_xlabel("Chi Squared Index")
#     axes[i].set_ylabel("Cumulative Probability")
#     axes[i].legend()
#     axes[i].grid(True)

# Adjust layout and show plot
# plt.tight_layout()
# plt.show()
