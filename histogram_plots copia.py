import matplotlib.pyplot as plt
from astropy.io.votable import parse

# Function to extract the fourth column from a VOTable file
def extract_fourth_column(file_path):
    votable = parse(file_path)
    table = votable.get_first_table().to_table()
    return table.columns[3]  # Fourth column (index 3)

x = "_sa04"

simbad_files = [
    "Simbad_variable_objects/variables_gSDSS"+x+".vot",
    "Simbad_variable_objects/variables_iSDSS"+x+".vot",
    "Simbad_variable_objects/variables_zSDSS"+x+".vot",
    "Simbad_variable_objects/variables_rSDSS"+x+".vot",
    "Simbad_variable_objects/variables_J0660"+x+".vot",
    "Simbad_variable_objects/variables_uJAVA"+x+".vot",
]

gaia_files = [
    "Gaia_nonvariable_objects/nonvariables_gSDSS"+x+".vot",
    "Gaia_nonvariable_objects/nonvariables_iSDSS"+x+".vot",
    "Gaia_nonvariable_objects/nonvariables_zSDSS"+x+".vot",
    "Gaia_nonvariable_objects/nonvariables_rSDSS"+x+".vot",
    "Gaia_nonvariable_objects/nonvariables_J0660"+x+".vot",
    "Gaia_nonvariable_objects/nonvariables_uJAVA"+x+".vot",
]

# simbad_files = [
#     "Simbad_variable_objects/variables_g_i"+x+".vot",
#     "Simbad_variable_objects/variables_g_J"+x+".vot",
#     "Simbad_variable_objects/variables_g_r"+x+".vot",
#     "Simbad_variable_objects/variables_g_z"+x+".vot",
#     "Simbad_variable_objects/variables_i_J"+x+".vot",
#     "Simbad_variable_objects/variables_i_r"+x+".vot",
#     "Simbad_variable_objects/variables_i_z"+x+".vot",
#     "Simbad_variable_objects/variables_J_r"+x+".vot",
#     "Simbad_variable_objects/variables_J_z"+x+".vot",
#     "Simbad_variable_objects/variables_r_z"+x+".vot",
# ]

# gaia_files = [
#     "Gaia_nonvariable_objects/nonvariables_g_i"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_g_J"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_g_r"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_g_z"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_i_J"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_i_r"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_i_z"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_J_r"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_J_z"+x+".vot",
#     "Gaia_nonvariable_objects/nonvariables_r_z"+x+".vot",
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

y = [0, 1.3] # range 

axes[0, 0].hist(data_x1, bins=150, alpha=0.6, color='blue', label="Variables", density=True, range=y)
axes[0, 0].hist(data_y1, bins=150, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
axes[0, 0].set_title("gSDSS")
axes[0, 0].set_xlim(0, 1.1)
# axes[0, 0].set_ylim(0, 20)  
axes[0, 0].legend()

axes[0, 1].hist(data_x2, bins=150, alpha=0.6, color='blue', label="Variables", density=True, range=y)
axes[0, 1].hist(data_y2, bins=150, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
axes[0, 1].set_title("iSDSS")
axes[0, 1].set_xlim(0, 1.1)
# axes[0, 1].set_ylim(0, 20)
axes[0, 1].legend()

axes[0, 2].hist(data_x3, bins=150, alpha=0.6, color='blue', label="Variables", density=True, range=y)
axes[0, 2].hist(data_y3, bins=150, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
axes[0, 2].set_title("zSDSS")
axes[0, 2].set_xlim(0, 1.1)
# axes[0, 2].set_ylim(0, 20)
axes[0, 2].legend()

axes[1, 0].hist(data_x4, bins=150, alpha=0.6, color='blue', label="Variables", density=True, range=y)
axes[1, 0].hist(data_y4, bins=150, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
axes[1, 0].set_title("rSDSS")
axes[1, 0].set_xlim(0, 1.1)
# axes[1, 0].set_ylim(0, 20)
axes[1, 0].legend()

axes[1, 1].hist(data_x5, bins=150, alpha=0.6, color='blue', label="Variables", density=True, range=y)
axes[1, 1].hist(data_y5, bins=150, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
axes[1, 1].set_title("J0660")
axes[1, 1].set_xlim(0, 1.1)
# axes[1, 1].set_ylim(0, 20)
axes[1, 1].legend()

axes[1, 2].hist(data_x6, bins=150, alpha=0.6, color='blue', label="Variables", density=True, range=y)
axes[1, 2].hist(data_y6, bins=150, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
axes[1, 2].set_title("uJAVA")
axes[1, 2].set_xlim(0, 1.1)
# axes[1, 2].set_ylim(0, 20)
axes[1, 2].legend()

# # 10 plots
# # Create subplots (2 rows, 3 columns)
# fig, axes = plt.subplots(2, 5, figsize=(15, 6))
# fig.suptitle("Normalised histogram comparing index values for varibale and non variable objects ", fontsize=16)

# y = [-2, 20] # range 

# # Plot manually for each subplot (adjust bins, limits, etc. as needed)
# axes[0, 0].hist(data_x1, bins=150, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[0, 0].hist(data_y1, bins=150, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[0, 0].set_title("gSDSS & iSDSS")
# axes[0, 0].set_xlim(-1, 5)
# # axes[0, 0].set_ylim(0, 20)  
# axes[0, 0].legend()

# axes[0, 1].hist(data_x2, bins=150, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[0, 1].hist(data_y2, bins=150, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[0, 1].set_title("gSDSS & J0660")
# axes[0, 1].set_xlim(-1, 5)
# # axes[0, 1].set_ylim(0, 20)
# axes[0, 1].legend()

# axes[0, 2].hist(data_x3, bins=150, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[0, 2].hist(data_y3, bins=150, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[0, 2].set_title("gSDSS & rSDSS")
# axes[0, 2].set_xlim(-1, 5)
# # axes[0, 2].set_ylim(0, 20)
# axes[0, 2].legend()

# axes[0, 3].hist(data_x4, bins=150, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[0, 3].hist(data_y4, bins=150, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[0, 3].set_title("gSDSS & zSDSS")
# axes[0, 3].set_xlim(-1, 5)
# # axes[1, 0].set_ylim(0, 20)
# axes[0, 3].legend()

# axes[0, 4].hist(data_x5, bins=150, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[0, 4].hist(data_y5, bins=150, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[0, 4].set_title("iSDSS & J0660")
# axes[0, 4].set_xlim(-1, 5)
# # axes[1, 1].set_ylim(0, 20)
# axes[0, 4].legend()

# axes[1, 0].hist(data_x6, bins=150, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[1, 0].hist(data_y6, bins=150, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[1, 0].set_title("iSDSS & rSDSS")
# axes[1, 0].set_xlim(-1, 5)
# # axes[1, 2].set_ylim(0, 20)
# axes[1, 0].legend()

# axes[1, 1].hist(data_x7, bins=150, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[1, 1].hist(data_y7, bins=150, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[1, 1].set_title("iSDSS & zSDSS")
# axes[1, 1].set_xlim(-1, 5)
# # axes[1, 2].set_ylim(0, 20)
# axes[1, 1].legend()

# axes[1, 2].hist(data_x8, bins=150, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[1, 2].hist(data_y8, bins=150, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[1, 2].set_title("J0660 & rSDSS")
# axes[1, 2].set_xlim(-1, 5)
# # axes[1, 2].set_ylim(0, 20)
# axes[1, 2].legend()

# axes[1, 3].hist(data_x9, bins=150, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[1, 3].hist(data_y9, bins=150, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[1, 3].set_title("J0660 & zSDSS")
# axes[1, 3].set_xlim(-1, 5)
# # axes[1, 2].set_ylim(0, 20)
# axes[1, 3].legend()

# axes[1, 4].hist(data_x10, bins=150, alpha=0.6, color='blue', label="Variables", density=True, range=y)
# axes[1, 4].hist(data_y10, bins=150, alpha=0.6, color='red', label="Non-variables", density=True, range=y)
# axes[1, 4].set_title("rSDSS & zSDSS")
# axes[1, 4].set_xlim(-1, 5)
# # axes[1, 2].set_ylim(0, 20)
# axes[1, 4].legend()


# Adjust layout
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('Figures/histogram'+x+'.png')
plt.show()