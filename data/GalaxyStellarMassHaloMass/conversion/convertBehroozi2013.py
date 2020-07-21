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

output_filename = "Behroozi2013.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

data = np.loadtxt("../raw/Behroozi_2013_z0p1.txt")
M_vir = (10 ** data[:, 0]) * unyt.Solar_Mass
M_star_ratio = (10 ** data[:, 1]) * unyt.dimensionless
M_star = M_vir * M_star_ratio

# Correct M_vir --> M_200 (fit from EAGLE cosmology - cosmology dependence is very weak)
M_200 = M_vir / 1.2

# Meta-data
comment = (
    "Fit obtained directly from Peter Behroozi's webpage. "
    "Stellar Masses: Chabrier IMF, BC03 SPS model, Blanton et al. dust model (i.e., kcorrect)."
    "Cosmology: Omega_m = 0.27, ns = 0.95, Omega_b = 0.046, sigma_8 = 0.82, h = 0.7."
    "No cosmology correction needed. "
    "Halo masses have been corrected from M_vir to M_200_cr."
    "Shows the relation betweeen stellar mass and halo mass."
)
citation = "Behroozi et al. (2013)"
bibcode = "2013ApJ...770...57B"
name = "Fit to the stellar mass - halo mass relation at z=0."
plot_as = "line"
redshift = 0.0
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_200, scatter=None, comoving=True, description="Halo Mass ($M_{200, {\rm crit}}$)"
)
processed.associate_y(
    M_star, scatter=None, comoving=True, description="Galaxy Stellar Mass"
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
