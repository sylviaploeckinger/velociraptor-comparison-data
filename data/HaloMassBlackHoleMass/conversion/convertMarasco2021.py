from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Marasco2021.txt"
delimiter = None
half_mass = 1
log_mass = 0

output_filename = "Marasco2021_Data.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

processed = ObservationalData()

# Read the data (only those columns we need here)
raw = np.loadtxt(input_filename, delimiter=delimiter,
                 usecols=(2, 3, 4, 5))

M_BH = 10 ** raw[:, 0] * unyt.Solar_Mass
M_BH_low = 10 ** (raw[:, 0] - raw[:, 1]) * unyt.Solar_Mass
M_BH_high = 10 ** (raw[:, 0] + raw[:, 1]) * unyt.Solar_Mass

M_halo = 10 ** raw[:, 2] * unyt.Solar_Mass
M_halo_low = 10 ** (raw[:, 2] - raw[:, 3]) * unyt.Solar_Mass
M_halo_high = 10 ** (raw[:, 2] + raw[:, 3]) * unyt.Solar_Mass

# Define the scatter as offset from the mean value
x_scatter = unyt.unyt_array((M_halo - M_halo_low, M_halo_high - M_halo))
y_scatter = unyt.unyt_array((M_BH - M_BH_low, M_BH_high - M_BH))

comment = (
    "Masses are provided h-free and cosmology-independent, so no "
    "h-correction made. "
    "Masses are (mostly) determined dynamically, with some stallar "
    "masses obtained from K-band luminosities with a fixed conversion "
    "factor. Halo masses are defined as M200_crit."
)
citation = "Marasco et al. (2021)"
bibcode = "2021arXiv210510508M"
name = "Halo Mass-Black Hole Mass"
plot_as = "points"
redshift = 0.0
h = cosmology.h

processed.associate_x(
    M_halo, scatter=x_scatter, comoving=False,
    description="Halo mass"
)
processed.associate_y(
    M_BH,
    scatter=y_scatter,
    comoving=False,
    description="Black hole mass",
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
