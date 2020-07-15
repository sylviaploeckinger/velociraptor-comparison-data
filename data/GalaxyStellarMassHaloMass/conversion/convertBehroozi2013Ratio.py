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

output_filename = "Behroozi2013Ratio.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

data = np.loadtxt("../raw/Behroozi_2013_z0p1.txt")
M_vir = (10 ** data[:, 0]) * unyt.Solar_Mass
M_star_ratio = (10 ** data[:, 1]) * unyt.dimensionless
M_star_ratio_p = (10 ** (data[:, 1] + data[:, 2])) * unyt.dimensionless
M_star_ratio_m = (10 ** (data[:, 1] - data[:, 3])) * unyt.dimensionless

# Define the scatter as offset from the mean value
y_scatter = unyt.unyt_array(
    (M_star_ratio - M_star_ratio_m, M_star_ratio_p - M_star_ratio_m)
)

# Meta-data
comment = (
    "Fit obtained directly from Peter Behroozi's webpage. "
    "Stellar Masses: Chabrier IMF, BC03 SPS model, Blanton et al. dust model (i.e., kcorrect)."
    "Cosmology: Omega_m = 0.27, ns = 0.95, Omega_b = 0.046, sigma_8 = 0.82, h = 0.7."
    "No cosmology correction needed. "
    "Shows the ratio betweeen stellar mass and halo mass."
)
citation = "Behroozi et al. (2013)"
bibcode = "2013ApJ...770...57B"
name = "Fit to the stellar mass - stellar halo mass relation at z=0."
plot_as = "line"
redshift = 0.0
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_vir, scatter=None, comoving=True, description="Halo Mass (M_200_cr)"
)
processed.associate_y(
    M_star_ratio,
    scatter=y_scatter,
    comoving=True,
    description="Galaxy Stellar Mass / Halo Mass (M_200_cr)",
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
