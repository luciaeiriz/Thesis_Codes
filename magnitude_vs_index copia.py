import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io.votable import parse
from astropy.table import Table

# Paths
data_folder = 'Prueba'
output_folder = 'Figures'

# Specify the specific file to process
target_file = 'CepOB3_J0660_NPh_aperphot_ststar_v4_part09_chisquared.vot'
file_path = os.path.join(data_folder, target_file)

# Read the VOTable
votable = parse(file_path)

# Access the first table
if votable.get_first_table() is not None:
    table = Table(votable.get_first_table().array)

    # Check if the required columns exist in the table
    if 'chi_squared' in table.colnames and 'mean_magnitude' in table.colnames:
        # Extract the data for plotting
        chi_squared = table['chi_squared']
        mean_magnitude = table['mean_magnitude']

        # Plotting
        plt.figure(figsize=(8, 6))
        plt.scatter(mean_magnitude, chi_squared,color='blue', s=10, alpha=0.7)
        plt.title('Mean Magnitude vs Chi Squared (J0660)')
        plt.xlabel('Chi Squared')
        plt.ylabel('Mean Magnitude')
        plt.grid(True)
        plt.tight_layout()

        # Save the plot
        output_file = os.path.join(output_folder, 'J0660_filtered.png')
        plt.savefig(output_file)

        # Show the plot
        plt.show()
    else:
        print("The required columns 'chi_squared' and 'mean_magnitude' are not found in the VOTable.")
else:
    print("No tables found in the VOTable.")