from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Cosmology
h_sim = cosmology.h

input_filename = "../raw/Kirby_2013_ascii.dat"
delimiter = "\t"

output_filename = "Kirby2013_Data.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data
raw = np.loadtxt(input_filename, delimiter=delimiter)
M_star = 10 ** raw[:, 0] * unyt.Solar_Mass
M_star_lo = 10 ** (raw[:, 0] - raw[:, 1]) * unyt.Solar_Mass
M_star_hi = 10 ** (raw[:, 0] + raw[:, 2]) * unyt.Solar_Mass

Z_star = 10 ** raw[:, 3] * unyt.dimensionless  # Z/Z_sun
Z_star_lo = 10 ** (raw[:, 3] - raw[:, 4]) * unyt.dimensionless  # Z/Z_sun
Z_star_hi = 10 ** (raw[:, 3] + raw[:, 5]) * unyt.dimensionless  # Z/Z_sun

# Define the scatter as offset from the mean value
x_scatter = unyt.unyt_array((M_star - M_star_lo, M_star_hi - M_star))
y_scatter = unyt.unyt_array((Z_star - Z_star_lo, Z_star_hi - Z_star))

# Meta-data
comment = (
    "Data obtained assuming a Kroupa IMF and 12+log10(Fe/H)_sun = 7.52. "
    "The metallicity is expressed in units of solar metallicity. "
)
citation = "Kirby et al. (2013)"
bibcode = "2013ApJ...779..102K"
name = "Stellar mass - Stellar metallicity relation"
plot_as = "points"
redshift = 0.0
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star, scatter=x_scatter, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    Z_star, scatter=y_scatter, comoving=True, description="Stellar metallicity"
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
