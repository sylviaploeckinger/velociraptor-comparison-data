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

input_filename = "../raw/Eckert2016.dat"

output_filename = "Eckert2016.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data
raw = np.loadtxt(input_filename)
M_500 = unyt.unyt_array((0.70 / h_sim) * 10 ** raw[:, 0], units="Msun")
fb_500 = unyt.unyt_array( raw[:, 1] * (0.70 / h_sim) ** (2.5)) / (Omega_b / Omega_m)

# Meta-data
comment = (
    "Based on observations obtained with XMM-Newton. "
    "Data was corrected for the simulation's cosmology."
)
citation = "Eckert et al. (2016)"
bibcode = "2016A&A...592A..12E"
name = "The XXL Survey. XIII. Baryon content of the bright cluster sample"
plot_as = "points"
redshift = 0.1
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_500, scatter=None, comoving=True, description="Halo mass (M_500)"
)
processed.associate_y(
    fb_500, scatter=None, comoving=True, description="Gas fraction (<R_500)"
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
