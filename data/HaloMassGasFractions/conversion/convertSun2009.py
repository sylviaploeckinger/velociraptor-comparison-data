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

input_filename = "../raw/Sun2009.dat"

output_filenames = "Sun2009.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data
raw = np.loadtxt(input_filename)
M_500 = unyt.unyt_array((10 ** 13) * (0.73 / h_sim) * raw[:, 1], units="Msun")
error_M_500_p = unyt.unyt_array((10 ** 13) * (0.73 / h_sim) * raw[:, 2], units="Msun")
error_M_500_m = unyt.unyt_array((10 ** 13) * (0.73 / h_sim) * raw[:, 3], units="Msun")
fb_500 = unyt.unyt_array((0.73 / h_sim) ** 1.5 * raw[:, 4], units="dimensionless")
error_fb_500_p = unyt.unyt_array(
    (0.73 / h_sim) ** 1.5 * raw[:, 5], units="dimensionless"
)
error_fb_500_m = unyt.unyt_array(
    (0.73 / h_sim) ** 1.5 * raw[:, 6], units="dimensionless"
)

# Define the scatter as offset from the mean value
x_scatter = unyt.unyt_array((error_M_500_m, error_M_500_p))
y_scatter = unyt.unyt_array((error_fb_500_m, error_fb_500_p))

# Meta-data
comment = (
    "Gas and total mass profiles for 23 low-redshift, relaxed groups spanning "
    "a temperature range 0.7-2.7 keV, derived from Chandra. "
    "Data was corrected for the simulation's cosmology."
)
citation = "Sun et al. (2009)"
bibcode = "2009ApJ...693.1142S"
name = "Halo mass - Gas fraction relation from 43 low-redshift Chandra-observed groups."
plot_as = "points"
redshift = 0.1
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_500, scatter=x_scatter, comoving=True, description="Halo mass (M_500)"
)
processed.associate_y(
    fb_500, scatter=y_scatter, comoving=True, description="Gas fraction (<R_500)"
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
