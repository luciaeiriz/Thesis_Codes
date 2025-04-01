from astropy.io.votable import parse

# Define colors for different OTYPE_S categories
otype_colors = {
    "EclBin": "green",
    "RSCVnV*": "green",
    "BYDraV*": "green",
    "LongPeriodV*": "blue",
    "LongPeriodV*_Candidate": "blue",
    "YSO": "red",
    "YSO_Candidate": "red",
    "Star": "yellow", 
    "EmLine*": "red",
    "SB*": "green",
    "Be*": "orange",
    "TTauri*": "red",
    "ClassicalCep": "blue",
    "**": "yellow",
    "Mira": "blue",
    "RRLyrae": "blue",
    "HighPM*": "magenta",
    "HighMassXBin": "green",
    "C*": "blue",
    "Planet": "brown",
    "C*_Candidate": "blue",
    "Cepheid": "blue",
    "OrionV*": "red", 
    "PulsV*": "yellow",
    "delSctV*": "yellow",
    "PlanetaryNeb": "blue",
    "RGB*": "blue",
    "Ae*": "red",
    "X": "white"
}

# Load the VOTable
votable = parse("Variable_objects/Stetson_Variables/matched_variables_stetson.vot")  # Update with actual file path
table = votable.get_first_table().to_table()

# Extract RA, DEC, and OTYPE_S columns
ra = table['RA']
dec = table['DEC']
otypes = table['OTYPE_S']

# Create DS9 region file
output_file = "ds9_regions.reg"

with open(output_file, "w") as f:
    f.write("# Region file format: DS9 version 4.1\n")
    f.write("fk5\n")  # Use fk5 for J2000 sky coordinates

    for r, d, otype in zip(ra, dec, otypes):
        color = otype_colors.get(otype, "red")  # Default color is red if type is unknown
        f.write(f"circle({r},{d},10\") # color={color}\n")  # 10" radius

print(f"DS9 region file created: {output_file}")






