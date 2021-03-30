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
M_HI = unyt.unyt_array([raw[2]] * 2, unyt.Solar_Mass)
M_H2 = unyt.unyt_array([raw[4]] * 2, unyt.Solar_Mass)
M_neut = M_HI + M_H2
MH2_p_Mneut = M_H2 / M_neut
Mneut_err = unyt.unyt_array(np.sqrt(raw[5] ** 2 + raw[3] ** 2), unyt.Solar_Mass)
MH2_err = unyt.unyt_array(raw[5], unyt.Solar_Mass)
MH2_p_Mneut_err = (
    np.sqrt((MH2_err / M_H2) ** 2 + (Mneut_err / M_neut) ** 2) * MH2_p_Mneut
)

x_scatter = unyt.unyt_array([raw[1] * unyt.Solar_Mass] * 2)
y_scatter = unyt.unyt_array([MH2_p_Mneut_err] * 2)

# Meta-data
comment = (
    "Stellar Masses obtained assuming a Chabrier (2003) IMF. "
    "local measurements decoupled from the Hubble flow (no h)."
    "Neutral measurements via combining NGVS CO detection, with"
    "HI data from the VLA archive."
)

citation = "VCC2062 (Lisenfeld et al. 2016)"
bibcode = "2016A&A...590A..92L"
name = "Stellar mass - H2 to Neutral Gas Mass ratio"
plot_as = "points"
redshift = 0.0
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star, scatter=x_scatter, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    MH2_p_Mneut,
    scatter=y_scatter,
    comoving=True,
    description="Stellar mass - H2 to Neutral Gas Mass ratio",
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
