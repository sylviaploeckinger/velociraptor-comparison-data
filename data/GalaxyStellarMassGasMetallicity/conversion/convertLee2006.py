from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Cosmologies
h_sim = cosmology.h

input_filename = "../raw/Lee_2006_ascii.txt"
delimiter = "\t"

output_filename = "Lee2006_Data.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data
raw = np.loadtxt(input_filename, delimiter=delimiter)
M_star = 10 ** raw[:, 0] * unyt.Solar_Mass
Z_median = raw[:, 1] * unyt.dimensionless  # 12 + log(O/H)
Z_sigma = raw[:, 2] * unyt.dimensionless  # 12 + log(O/H)

# Define the scatter as offset from the mean value
y_scatter = unyt.unyt_array((Z_sigma, Z_sigma))

# Convert the metallicities to the Tremonti+2004 convention
Z_median = Z_median + (8.69 - 8.66)

# Meta-data
comment = (
    "Data obtained assuming a Chabrier IMF and 12+log10(O/H)_sun = 8.66. "
    "The metallicity is expressed as 12 + log10(O/H). "
    "In these units the solar metallicity is 8.69."
    "There is also a correction to convert the metallicities to the 12+log10(O/H) convention."
)
citation = "Lee et al. (2006)"
bibcode = "2006ApJ...647..970L"
name = "Stellar mass - Gas phase metallicity relation"
plot_as = "points"
redshift = 0.0
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star, scatter=None, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    Z_median, scatter=y_scatter, comoving=True, description="Gas phase metallicity"
)
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(redshift)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)
