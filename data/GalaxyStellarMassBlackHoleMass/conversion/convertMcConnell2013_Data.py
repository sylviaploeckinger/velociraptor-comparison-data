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

input_filename = "../raw/McConnell_Ma_2013_ascii.txt"
delimiter = "\t"

output_filename = "McConnell2013_Data.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data
raw = np.loadtxt(input_filename, delimiter=delimiter, usecols=(2, 3, 4, 11))
M_BH = raw[:, 0] * unyt.Solar_Mass
M_BH_lo = raw[:, 1] * unyt.Solar_Mass
M_BH_up = raw[:, 2] * unyt.Solar_Mass
M_star = raw[:, 3] * unyt.Solar_Mass

# Only keep the data for which the stellar mass was measured
mask = M_star > 0
M_BH = M_BH[mask]
M_BH_lo = M_BH_lo[mask]
M_BH_up = M_BH_up[mask]
M_star = M_star[mask]

# Apply the 0.24 dex undertainty on the stellar mass
M_star_up = (10.0 ** (np.log10(M_star / unyt.Solar_Mass) + 0.24)) * unyt.Solar_Mass
M_star_lo = (10.0 ** (np.log10(M_star / unyt.Solar_Mass) - 0.24)) * unyt.Solar_Mass

# Define the scatter as offset from the mean value
x_scatter = unyt.unyt_array((M_star - M_star_lo, M_star_up - M_star))
y_scatter = unyt.unyt_array((M_BH - M_BH_lo, M_BH_up - M_BH))

# Meta-data
comment = (
    "Data obtained assuming the total stellar mass is the same as the bulge mass. "
    "No cosmology correction needed."
)
citation = "McConnell & Ma (2013) (Data)"
bibcode = "2013ApJ...764..184M"
name = "Black hole mass - stellar mass relation from 36 local galaxies."
plot_as = "points"
redshift = 0.0
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star, scatter=x_scatter, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    M_BH, scatter=y_scatter, comoving=True, description="Black Hole Mass"
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
