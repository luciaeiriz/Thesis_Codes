import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io.votable import parse
from astropy.table import Table

# Paths
data_folder_1 = 'Gaia_nonvariable_objects'
data_folder_2 = 'Simbad_variable_objects'
output_folder = 'Figures'

# Specify the files to process
index_name = 'chisquared' #change
filter = 'uJAVA'
file_1 = 'nonvariables_'+filter+'_'+index_name+'.vot'
file_2 = 'variables_'+filter+'_'+index_name+'.vot'

# Read the first VOTable (Non-variable objects)
file_path_1 = os.path.join(data_folder_1, file_1)
votable_1 = parse(file_path_1)
table_1 = Table(votable_1.get_first_table().array)

# Read the second VOTable (Variable objects)
file_path_2 = os.path.join(data_folder_2, file_2)
votable_2 = parse(file_path_2)
table_2 = Table(votable_2.get_first_table().array)

# Check if required columns exist #change
if 'chi_squared' in table_1.colnames and 'mean_magnitude' in table_1.colnames and \
   'chi_squared' in table_2.colnames and 'mean_magnitude' in table_2.colnames:

    # Extract data
    index_1 = table_1['chi_squared'] #change    
    mean_magnitude_1 = table_1['mean_magnitude']
    
    index_2 = table_2['chi_squared'] #change  
    mean_magnitude_2 = table_2['mean_magnitude']

    # Plotting
    # plt.figure(figsize=(8, 6))
    # plt.scatter(mean_magnitude_1, index_1, color='blue', s=10, alpha=0.7, label="Non-variable objects")
    # plt.scatter(mean_magnitude_2, index_2, color='red', s=10, alpha=0.7, label="Variable objects")

    # plt.title('Mean Magnitude vs IQR (rSDSS)') #change 
    # plt.xlabel('Mean Magnitude')
    # plt.ylabel('Interquartile Range') #change 
    # plt.legend()
    # plt.grid(True)
    # plt.tight_layout()

    plt.hist(mean_magnitude_1, bins=70,density=True ,color='red', alpha=0.6, label="Non-variable objects")
    plt.hist(mean_magnitude_2, bins=70,density=True ,color='blue', alpha=0.5, label="Variable objects")

    plt.title(filter)
    plt.xlabel('Mean Magnitude', fontsize=14)
    plt.ylabel('Count', fontsize=14)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    # Save the plot
    output_file = os.path.join(output_folder, 'mean_vs_index_'+index_name+'.png')
    # plt.savefig(output_file)

    plt.savefig('Figures/mean_magnitude_'+ filter+'.png')
else:
    print("The required columns 'median_absolute_deviation' and 'mean_magnitude' are not found in one or both VOTables.")


    