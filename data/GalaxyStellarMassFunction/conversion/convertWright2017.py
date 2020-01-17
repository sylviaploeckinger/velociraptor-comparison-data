from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Cosmologies
h_obs = 0.7
h_sim = cosmology.h

input_filename = "../raw/Wright2017.txt"
delimiter = " "

output_filename = "Wright2017.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

processed = ObservationalData()
raw = np.loadtxt(input_filename, delimiter=delimiter)

log10_Mstar = raw[:, 0]
Phi = raw[:, 1]
Phi_minus = raw[:, 2]
Phi_plus = raw[:, 3]

# Convert to SWIFT-friendly units and correct for cosmology
Mstar = (10.0 ** log10_Mstar) * unyt.Solar_Mass * (h_sim / h_obs) ** -2
Phi = Phi * unyt.Mpc ** (-3) * (h_sim / h_obs) ** 3
Phi_plus = (Phi_plus * unyt.Mpc ** (-3)) * (h_sim / h_obs) ** 3
Phi_minus = (Phi_minus * unyt.Mpc ** (-3)) * (h_sim / h_obs) ** 3

# Meta-data
comment = f"Data obtained assuming a Chabrier IMF and h = 0.7. h-corrected for SWIFT using cosmology: {cosmology.name}. Ignoring the mass bins for which GAMA is systematically incomplete."
citation = "Wright et al. (2017) (GAMA-II)"
bibcode = "2017MNRAS.470.283W"
name = "GSMF from GAMA-II using bivariate brightness distribution"
redshift = 0.1
plot_as = "points"

# Write everything
processed.associate_x(
    Mstar, scatter=None, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    Phi,
    scatter=unyt.unyt_array([Phi_minus, Phi_plus]),
    comoving=True,
    description="Phi (GSMF)",
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
