import numpy as np
import matplotlib.pyplot as plt
from astropy.io.votable import parse
import os

# Function to plot histogram of the fourth column with values below a certain threshold
def plot_histogram(vot_file):
    # Read the VOTable
    vot_data = parse(vot_file)
    
    # Assuming the data is in the first table of the VOTable
    table = vot_data.get_first_table().array
    
    # Extract the fourth column values (chi_squared)
    chi_squared = table[table.dtype.names[3]]  # Accessing the fourth column by index

    # Filter chi_squared to only include values below 250
    chi_squared_filtered = chi_squared[chi_squared < 10]

    # Create a histogram of the filtered 'chi_squared' values
    plt.figure(figsize=(10, 6))
    plt.hist(chi_squared_filtered, bins=50, color='blue', alpha=0.7, edgecolor='black')
    
    plt.title('Histogram of Chi-squared Values (Filtered)')
    plt.xlabel('Chi-squared')
    plt.ylabel('Counts in Bin')
    
    # # Set the y-axis to logarithmic scale
    # plt.yscale('log')

    # Add grid and other formatting
    plt.grid()
    
    # Ensure the 'Figures' folder exists
    if not os.path.exists('Figures'):
        os.makedirs('Figures')
    
    # Save the figure
    plt.savefig('Figures/gSDSS_part01_chisqaured_pernight_filtered.png')
    plt.close()  # Close the figure to free up memory

# Specify the path to your VOTable file
vot_file_path = 'Chi_squared_Indices/CepOB3_gSDSS_NPh_aperphot_ststar_v4_part01_chisquared.vot'

# Call the function to plot the histogram
plot_histogram(vot_file_path)

print("Filtered histogram saved in the 'Figures' folder.")
