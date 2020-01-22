from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Baldry2012.txt"
delimiter = None

output_filename = "Baldry2012.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

processed = ObservationalData()
raw = np.loadtxt(input_filename, delimiter=delimiter)

comment = (
    "Baldry et al 2012, Table 1, Chabrier IMF. Below 10^8 solar mass, these are "
    "lower limits. Total stellar mass density of 2.3e8 Msun / Mpc^3. Poisson errors. "
    f"h-corrected for SWIFT using Cosmology: {cosmology.name}."
)
citation = "Baldry et al. (2012)"
bibcode = "2012MNRAS.421..621B"
name = "GSMF from GAMA"
plot_as = "points"
redshift = 0.06
h = cosmology.h
h_obs = 0.7

log_M = raw.T[0]
M = 10 ** (log_M) * unyt.Solar_Mass * (h / h_obs) ** -2
Phi = raw.T[2] * unyt.Mpc ** (-3) * (h / h_obs) ** 3
Phi_scatter = raw.T[3] * unyt.Mpc ** (-3) * (h / h_obs) ** 3

processed.associate_x(M, scatter=None, comoving=True, description="Galaxy Stellar Mass")
processed.associate_y(Phi, scatter=Phi_scatter, comoving=True, description="Phi, GSMF")
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
