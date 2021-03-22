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

input_filename = "../raw/Grossi2016.txt"

output_filename = "Grossi2016_Data.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data
raw = np.loadtxt(input_filename)
M_star = pow(10.0, raw[:, 1]) * unyt.Solar_Mass
M_H2 = pow(10.0, raw[:, 2]) * unyt.Solar_Mass
M_H2_perr = pow(10.0, raw[:, 2] + raw[:, 3]) * unyt.Solar_Mass
M_H2_merr = pow(10.0, raw[:, 2] - raw[:, 3]) * unyt.Solar_Mass

# index galaxies with H2 detections
index_detections = raw[:, 2] > -99
M_star = M_star[index_detections]
M_H2_frac = M_H2[index_detections] / M_star
M_H2_frac_perr = M_H2_perr[index_detections] / M_star
M_H2_frac_merr = M_H2_merr[index_detections] / M_star

y_scatter = unyt.unyt_array([M_H2_frac - M_H2_frac_merr, M_H2_frac_perr - M_H2_frac])

# Meta-data
comment = (
    "Stellar Masses obtained assuming a Kroupa (2001) IMF. "
    "local measurements decoupled from the Hubble flow (no h)."
    "H2 measurements via CO detections in the Virgo Cluster "
    "Catalogue."
)

citation = "Grossi et al 2016 (CO, VCC)"
bibcode = "2016A&A...590A..27G"
name = "Stellar mass - H2 Gas to Stellar Mass fractions"
plot_as = "points"
redshift = 0.0
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star, scatter=None, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    M_H2_frac,
    scatter=y_scatter,
    comoving=True,
    description="H2 Gas to Stellar Mass Fraction",
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
