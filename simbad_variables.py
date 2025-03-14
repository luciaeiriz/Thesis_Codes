from astropy.table import Table
import os

# Define the input and output file paths
input_file = 'simbad.vot'
output_file = 'simbad_variables.vot'

# List of variable star types to filter
variable_types = [
    "LongPeriodV*", "LongPeriodV*_Candidate", "EclBin", "RSCVnV*", "BYDraV*", "RRLyrae", "bCepV*", "ClassicalCep", "Mira", 
    "PulsV*", "Cepheid", "delSctV*", "Variable*", "OrionV*", "CataclyV*", "SB*", "YSO", "LongPeriodV*_Candidate", "YSO_Candidate",
    "EmLine*", "HerbigHaroObj", "TTauri*", "RRLyrae_Candidate", "Cepheid", "BrownD*_Candidate"
]

# Load the VOTable
try:
    table = Table.read(input_file, format='votable')
except Exception as e:
    print(f"Error reading the VOTable: {e}")
    exit()

# Extract rows where the first column contains any of the variable types
try:
    filtered_table = table[[row[0] in variable_types for row in table]]
except KeyError:
    print("Error: Unable to filter rows. Check if the first column contains valid data.")
    exit()

# Save the filtered table to a new VOTable
try:
    filtered_table.write(output_file, format='votable', overwrite=True)
    print(f"Filtered VOTable saved to {output_file}")
except Exception as e:
    print(f"Error writing the filtered VOTable: {e}")
