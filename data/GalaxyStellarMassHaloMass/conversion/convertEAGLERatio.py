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

output_filename = "Schaye2015_Ref_100_z0p0_Ratio.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

data = np.loadtxt("../raw/EAGLE_SMHM_ratio_L100N1504_z0p0.txt")
M_200 = data[:, 0] * unyt.Solar_Mass
M_star_ratio_16 = data[:, 1] * unyt.dimensionless
M_star_ratio_50 = data[:, 2] * unyt.dimensionless
M_star_ratio_84 = data[:, 3] * unyt.dimensionless

# Define the scatter as offset from the mean value
y_scatter = unyt.unyt_array(
    (M_star_ratio_50 - M_star_ratio_16, M_star_ratio_84 - M_star_ratio_50)
)

# Meta-data
comment = (
    "Medians obtained directly from the subfind catalogs. "
    "Central galaxies only. Halos with 0 stellar mass included. "
    "No cosmology correction needed."
)
citation = "EAGLE - L100N1504"
bibcode = "2015MNRAS.446..521S"
name = "Galaxy Stellar Mass / Halo mass - Halo Mass relation at z=0."
plot_as = "line"
redshift = 0.0
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_200, scatter=None, comoving=True, description="Halo Mass ($M_{200, {\rm crit}}$)"
)
processed.associate_y(
    M_star_ratio_50,
    scatter=y_scatter,
    comoving=True,
    description="Galaxy Stellar Mass (30kpc, 3D) / Halo mass ($M_{200, {\rm crit}}$)",
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
