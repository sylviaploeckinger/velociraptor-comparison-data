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
Omega_b = cosmology.Ob0
Omega_m = cosmology.Om0

input_filename = "../raw/Lin2012.dat"

output_filename = "Lin2012.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data
raw = np.loadtxt(input_filename)
M_500 = unyt.unyt_array((0.71 / h_sim) * 10 ** raw[:, 0], units="Msun")
M_500_error = unyt.unyt_array((0.71 / h_sim) * raw[:, 1], units="Msun")
M_500_gas = unyt.unyt_array((0.71 / h_sim) * 10 ** raw[:, 2], units="Msun")
M_500_gas_error = unyt.unyt_array((0.71 / h_sim) * raw[:, 3], units="Msun")
z = raw[:, 6]

# Compute the gas fractions
fb_500 = (M_500_gas / M_500) * (0.71 / h_sim) ** (2.5)
fb_500_error = fb_500 * ((M_500_error / M_500) + (M_500_gas_error / M_500_gas))

# Normalise by the cosmic mean
fb_500 = fb_500 / (Omega_b / Omega_m)
fb_500_error = fb_500_error / (Omega_b / Omega_m)

# Select only the low-z data
M_500 = M_500[z < 0.25]
fb_500 = fb_500[z < 0.25]
M_500_error = M_500_error[z < 0.25]
fb_500_error = fb_500_error[z < 0.25]

# Define the scatter as offset from the mean value
x_scatter = M_500_error
y_scatter = fb_500_error

# Meta-data
comment = (
    "Ionized gas out of 94 clusters combining Chandra, WISE and 2MASS. "
    "Data was corrected for the simulation's cosmology."
)
citation = "Lin et al. (2012) - (z < 0.25)"
bibcode = "2012ApJ...745L...3L"
name = "Halo mass - Gas fraction relation from Chandra-observed clusters."
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
