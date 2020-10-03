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

input_filename = "../raw/Vikhlinin2006.dat"

output_filenames = "Vikhlinin2006.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data
raw = np.loadtxt(input_filename)
M_500 = unyt.unyt_array((10 ** 14) * (0.72 / h_sim) * raw[:, 1], units="Msun")
error_M_500 = unyt.unyt_array((10 ** 14) * (0.72 / h_sim) * raw[:, 2], units="Msun")
fb_500 = unyt.unyt_array((0.72 / h_sim) ** 1.5 * raw[:, 3], units="dimensionless")
error_fb_500 = unyt.unyt_array((0.72 / h_sim) ** 1.5 * raw[:, 4], units="dimensionless")

# Define the scatter as offset from the mean value
x_scatter = unyt.unyt_array((error_M_500, error_M_500))
y_scatter = unyt.unyt_array((error_fb_500, error_fb_500))

# Meta-data
comment = (
    "Gas and total mass profiles for 13 low-redshift, relaxed clusters spanning "
    "a temperature range 0.7-9 keV, derived from all available Chandra data of "
    "sufficient quality. Data was corrected for the simulation's cosmology."
)
citation = "Vikhlinin et al. (2006)"
bibcode = "2006ApJ...640..691V"
name = "Halo mass - Gas fraction relation from 13 low-redshift Chandra-observed relaxed clusters."
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
