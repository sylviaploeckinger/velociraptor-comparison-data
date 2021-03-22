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

input_filename = "../raw/Gratier2010.txt"

output_filename = "Gratier2010_Data.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data
# No error as error is only stated for H2, making any propagation an underestimate.
raw = np.loadtxt(input_filename)
M_star = unyt.unyt_array([raw[0]] * 2, unyt.Solar_Mass)
M_HI = unyt.unyt_array([raw[2]] * 2, unyt.Solar_Mass)
M_H2 = unyt.unyt_array([raw[4]] * 2, unyt.Solar_Mass)
MH2_p_Mneut = M_H2 / (M_H2 + M_HI)

# Meta-data
comment = "Local measurements decoupled from the Hubble flow (no h)."

citation = "NGC6822 (Gratier et al. 2010)"
bibcode = "2010A&A...512A..68G"
name = "Stellar mass - H2 to Neutral Gas ratio"
plot_as = "points"
redshift = 0.0
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star, scatter=None, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    MH2_p_Mneut,
    scatter=None,
    comoving=True,
    description="Stellar mass - H2 to Neutral Gas ratio",
)
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(redshift, 0, 2)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)
