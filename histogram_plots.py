import matplotlib.pyplot as plt
from astropy.io.votable import parse
import matplotlib.ticker as ticker


# Function to extract the fourth column from a VOTable file
def extract_fourth_column(file_path):
    votable = parse(file_path)
    table = votable.get_first_table().to_table()
    return table.columns[3]  # Fourth column (index 3)

x = "l1"

simbad_files = [
    "Simbad_variable_objects/variables_gSDSS_"+x+".vot",
    "Simbad_variable_objects/variables_iSDSS_"+x+".vot",
    "Simbad_variable_objects/variables_zSDSS_"+x+".vot",
    "Simbad_variable_objects/variables_rSDSS_"+x+".vot",
    "Simbad_variable_objects/variables_J0660_"+x+".vot",
    "Simbad_variable_objects/variables_uJAVA_"+x+".vot",
]

gaia_files = [
    "Gaia_nonvariable_objects/nonvariables_gSDSS_"+x+".vot",
    "Gaia_nonvariable_objects/nonvariables_iSDSS_"+x+".vot",
    "Gaia_nonvariable_objects/nonvariables_zSDSS_"+x+".vot",
    "Gaia_nonvariable_objects/nonvariables_rSDSS_"+x+".vot",
    "Gaia_nonvariable_objects/nonvariables_J0660_"+x+".vot",
    "Gaia_nonvariable_objects/nonvariables_uJAVA_"+x+".vot",
]

# simbad_files = [
#     "Simbad_variable_objects/variables_g_i_"+x+".vot",
#     "Simbad_variable_objects/variables_g_J_"+x+".vot",
#     "Simbad_variable_objects/variables_g_r_"+x+".vot",
#     "Simbad_variable_objects/variables_g_z_"+x+".vot",
#     "Simbad_variable_objects/variables_i_J_"+x+".vot",
#     "Simbad_variable_objects/variables_i_r_"+x+".vot",
#     "Simbad_variable_objects/variables_i_z_"+x+".vot",
#     "Simbad_variable_objects/variables_J_r_"+x+".vot",
#     "Simbad_variable_objects/variables_J_z_"+x+".vot",
#     "Simbad_variable_objects/variables_r_z_"+x+".vot",
# ]

# gaia_files = [
#     "Gaia_nonvariable_objects/nonvariables_g_i_"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_g_J_"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_g_r_"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_g_z_"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_i_J_"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_i_r_"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_i_z_"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_J_r_"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_J_z_"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_r_z_"+x+".vot",
# ]

# simbad_files = [
#     "Simbad_variable_objects/variables_g_i_"+x+".vot",
#     "Simbad_variable_objects/variables_g_r_"+x+".vot",
#     "Simbad_variable_objects/variables_g_z_"+x+".vot",
#     "Simbad_variable_objects/variables_r_i_"+x+".vot",
#     "Simbad_variable_objects/variables_i_z_"+x+".vot",
#     "Simbad_variable_objects/variables_r_z_"+x+".vot",
# ]

# gaia_files = [
#     "Gaia_nonvariable_objects/nonvariables_g_i_"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_g_r_"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_g_z_"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_r_i_"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_i_z_"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_r_z_"+x+".vot",
# ]

data_x1 = extract_fourth_column(simbad_files[0])
data_y1 = extract_fourth_column(gaia_files[0])
data_x2 = extract_fourth_column(simbad_files[1])
data_y2 = extract_fourth_column(gaia_files[1])
data_x3 = extract_fourth_column(simbad_files[2])
data_y3 = extract_fourth_column(gaia_files[2])
data_x4 = extract_fourth_column(simbad_files[3])
data_y4 = extract_fourth_column(gaia_files[3])
data_x5 = extract_fourth_column(simbad_files[4])
data_y5 = extract_fourth_column(gaia_files[4])
data_x6 = extract_fourth_column(simbad_files[5])
data_y6 = extract_fourth_column(gaia_files[5])
# data_x7 = extract_fourth_column(simbad_files[6])
# data_y7 = extract_fourth_column(gaia_files[6])
# data_x8 = extract_fourth_column(simbad_files[7])
# data_y8 = extract_fourth_column(gaia_files[7])
# data_x9 = extract_fourth_column(simbad_files[8])
# data_y9 = extract_fourth_column(gaia_files[8])
# data_x10 = extract_fourth_column(simbad_files[9])
# data_y10 = extract_fourth_column(gaia_files[9])


# 6 plots
# Create subplots (2 rows, 3 columns)
fig, axes = plt.subplots(2, 3, figsize=(12, 8))
fig.suptitle("Normalised histogram comparing index values for varibale and non variable objects ", fontsize=16)

y = [-1, 1] # range 
index = "Lag-1"

axes[0, 0].hist(data_x1, bins=75, alpha=0.6, color='blue', label="Variables", density=True, range=y)
axes[0, 0].hist(data_y1, bins=75, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
axes[0, 0].set_title("gSDSS")
axes[0, 0].set_xlim(-0.5, 1)
#axes[0, 0].set_ylim(0.25, 2)  
axes[0, 0].legend()

axes[0, 1].hist(data_x2, bins=75, alpha=0.6, color='blue', label="Variables", density=True, range=y)
axes[0, 1].hist(data_y2, bins=75, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
axes[0, 1].set_title("gSDSS & rSDSS")
#axes[0, 1].set_xlim(0, 0.04)
#axes[0, 1].set_ylim(0, 2)
axes[0, 1].legend()

axes[0, 2].hist(data_x3, bins=75, alpha=0.6, color='blue', label="Variables", density=True, range=y)
axes[0, 2].hist(data_y3, bins=75, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
axes[0, 2].set_title("gSDSS & zSDSS")
#axes[0, 2].set_xlim(0, 0.03)
#axes[0, 2].set_ylim(0, 1.4)
axes[0, 2].legend()

axes[1, 0].hist(data_x4, bins=75, alpha=0.6, color='blue', label="Variables", density=True, range=y)
axes[1, 0].hist(data_y4, bins=75, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
axes[1, 0].set_title("iSDSS & rSDSS")
#axes[1, 0].set_xlim(0, 0.04)
#axes[1, 0].set_ylim(0, 0.7)
axes[1, 0].legend()

axes[1, 1].hist(data_x5, bins=75, alpha=0.6, color='blue', label="Variables", density=True, range=y)
axes[1, 1].hist(data_y5, bins=75, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
axes[1, 1].set_title("J0660")
axes[1, 1].set_xlim(-0.5, 1)
#axes[1, 1].set_ylim(0, 0.6)
axes[1, 1].legend()
axes[1, 1].set_xlabel(index)

axes[1, 2].hist(data_x6, bins=75, alpha=0.6, color='blue', label="Variables", density=True, range=y)
axes[1, 2].hist(data_y6, bins=75, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
axes[1, 2].set_title("rSDSS & zSDSS")
#axes[1, 2].set_xlim(0, 0.024)
#axes[1, 2].set_ylim(0, 0.4)
axes[1, 2].legend()

# # 10 plots
# # Create subplots (2 rows, 3 columns)
# fig, axes = plt.subplots(2, 5, figsize=(15, 6))
# fig.suptitle("Normalised histogram comparing index values for varibale and non variable objects ", fontsize=16)

# y = [-0.5, 10] # range 

# # Plot manually for each subplot (adjust bins, limits, etc. as needed)
# axes[0, 0].hist(data_x1, bins=100, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[0, 0].hist(data_y1, bins=100, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[0, 0].set_title("gSDSS & iSDSS")
# axes[0, 0].set_xlim(-0.5, 4)
# axes[0, 0].set_ylim(0, 4)  
# axes[0, 0].legend()

# axes[0, 1].hist(data_x2, bins=100, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[0, 1].hist(data_y2, bins=100, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[0, 1].set_title("gSDSS & J0660")
# axes[0, 1].set_xlim(-0.5, 4)
# axes[0, 1].set_ylim(0, 4)
# axes[0, 1].legend()

# axes[0, 2].hist(data_x3, bins=100, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[0, 2].hist(data_y3, bins=100, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[0, 2].set_title("gSDSS & rSDSS")
# axes[0, 2].set_xlim(-0.5, 4)
# axes[0, 2].set_ylim(0, 4)
# axes[0, 2].legend()

# axes[0, 3].hist(data_x4, bins=100, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[0, 3].hist(data_y4, bins=100, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[0, 3].set_title("gSDSS & zSDSS")
# axes[0, 3].set_xlim(-0.5, 4)
# axes[1, 0].set_ylim(0, 4)
# axes[0, 3].legend()

# axes[0, 4].hist(data_x5, bins=100, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[0, 4].hist(data_y5, bins=100, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[0, 4].set_title("iSDSS & J0660")
# axes[0, 4].set_xlim(-0.5, 4)
# axes[1, 1].set_ylim(0, 4)
# axes[0, 4].legend()

# axes[1, 0].hist(data_x6, bins=100, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[1, 0].hist(data_y6, bins=100, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[1, 0].set_title("iSDSS & rSDSS")
# axes[1, 0].set_xlim(-0.5, 4)
# axes[1, 2].set_ylim(0, 4)
# axes[1, 0].legend()

# axes[1, 1].hist(data_x7, bins=100, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[1, 1].hist(data_y7, bins=100, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[1, 1].set_title("iSDSS & zSDSS")
# axes[1, 1].set_xlim(-0.5, 4)
# axes[1, 2].set_ylim(0, 4)
# axes[1, 1].legend()

# axes[1, 2].hist(data_x8, bins=100, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[1, 2].hist(data_y8, bins=100, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[1, 2].set_title("J0660 & rSDSS")
# axes[1, 2].set_xlim(-0.5, 4)
# axes[1, 2].set_ylim(0, 4)
# axes[1, 2].legend()

# axes[1, 3].hist(data_x9, bins=100, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[1, 3].hist(data_y9, bins=100, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[1, 3].set_title("J0660 & zSDSS")
# axes[1, 3].set_xlim(-0.5, 4)
# axes[1, 2].set_ylim(0, 4)
# axes[1, 3].legend()

# axes[1, 4].hist(data_x10, bins=100, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[1, 4].hist(data_y10, bins=100, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[1, 4].set_title("rSDSS & zSDSS")
# axes[1, 4].set_xlim(-0.5, 4)
# axes[1, 2].set_ylim(0, 4)
# axes[1, 4].legend()


# Adjust layout
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()