from astropy.io.votable import parse

# Load the VOTable
votable = parse("Variable_objects/Chi_Squared_Variables/variables_chisquared.vot")  # Update with your actual file path
table = votable.get_first_table().to_table()

# Extract RA & DEC columns (already in degrees)
ra = table['RA']
dec = table['DEC']

# Create DS9 region file
output_file = "ds9_regions.reg"

with open(output_file, "w") as f:
    f.write("# Region file format: DS9 version 4.1\n")
    f.write("fk5\n")  # Use fk5 for J2000 sky coordinates

    for r, d in zip(ra, dec):
        f.write(f"circle({r},{d},10\") # color=red\n")  # 10" radius, red color

print(f"DS9 region file created: {output_file}")





