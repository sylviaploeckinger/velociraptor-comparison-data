from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import re
import sys
import itertools as it

ORIGINAL_H = 0.7

# Exec the master cosmology file passed as first argument
# These lines are _required_ and you are required to use
# the cosmology specified (this is an astropy.cosmology
# instance)
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Jones2018.txt"
output_filename = "Jones2018.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

data = np.genfromtxt(input_filename, comments="#")

processed = ObservationalData()

comment = (
    "Jones et al. (2018). Estimated for ALFALFA survey at z=0"
    f"data, h-corrected for SWIFT using Cosmology: {cosmology.name}."
)

citation = "Jones et al. (2018), z = 0"
bibcode = "2018MNRAS.477....2J"
name = "HIMF from ALFALFA at z=0"
plot_as = "points"
redshift = 0.0
h = cosmology.h

M = 10 ** data[:, 0] * (h / ORIGINAL_H) ** (-2) * unyt.Solar_Mass
Phi = 10 ** data[:, 1] * (h / ORIGINAL_H) ** 3 * unyt.Mpc ** (-3)

# no error in M_HI provided
M_err = np.row_stack([M * 0.0] * 2) * M.units
Phi_err = np.abs(
    (10 ** data[:, [2, 3]] * (h / ORIGINAL_H) ** 3 * unyt.Mpc ** (-3)) - Phi[:, None]
).T

processed.associate_x(M, scatter=M_err, comoving=True, description="Galaxy HI Mass")
processed.associate_y(Phi, scatter=Phi_err, comoving=True, description="Phi (HIMF)")
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
