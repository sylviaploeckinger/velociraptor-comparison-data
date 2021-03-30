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

input_filename = "../raw/Lisenfeld2016.txt"

output_filename = "Lisenfeld2016_Data.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data
raw = np.loadtxt(input_filename)
M_star = unyt.unyt_array([raw[0]] * 2, unyt.Solar_Mass)
M_H2 = unyt.unyt_array([raw[4]] * 2, unyt.Solar_Mass)
MH2_p_Mstar = M_H2 / M_star
MH2_p_Mstar_err = np.sqrt((raw[5] / raw[4]) ** 2 + (raw[1] / raw[0]) ** 2) * MH2_p_Mstar

x_scatter = unyt.unyt_array([raw[1] * unyt.Solar_Mass] * 2)
y_scatter = unyt.unyt_array([MH2_p_Mstar_err] * 2)

# Meta-data
comment = (
    "Stellar Masses obtained assuming a Chabrier (2003) IMF. "
    "local measurements decoupled from the Hubble flow (no h)."
    "H2 measurements via CO detection from NGVS."
)

citation = "VCC2062 (Lisenfeld et al. 2016)"
bibcode = "2016A&A...590A..92L"
name = "Stellar mass - H2 Gas to Stellar Mass ratio"
plot_as = "points"
redshift = 0.0
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star, scatter=x_scatter, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    MH2_p_Mstar,
    scatter=y_scatter,
    comoving=True,
    description="Stellar mass - H2 Gas to Stellar Mass ratio",
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
