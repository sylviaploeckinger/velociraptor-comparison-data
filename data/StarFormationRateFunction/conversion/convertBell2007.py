from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

output_filename = "Bell2007.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# The data is given in the form of a Schechter fit (Table 1 of the paper)
phi_star = 9.8 * (10.0 ** -4) * unyt.Mpc ** (-3)  # Mpc^-3 dex^-1
alpha = -1.45
SFR_star = (10 ** 1.11) * unyt.Solar_Mass / unyt.year  # Msun  / year

# Minimal and maximal mass
SFR_min = (10 ** -1) * unyt.Solar_Mass / unyt.year  # Msun / year
SFR_max = (10 ** 2) * unyt.Solar_Mass / unyt.year  # Msun / year

# Create the x-data
SFR = (
    np.logspace(np.log10(SFR_min), np.log10(SFR_max), 50) * unyt.Solar_Mass / unyt.year
)  # Msun / year

# Create the y-data (Schechter function)
Phi = (phi_star / SFR_star) * np.exp(-SFR / SFR_star) * (SFR / SFR_star) ** alpha
Phi = Phi * SFR * np.log(10)

# Meta-data
comment = "Data obtained assuming a Chabrier IMF and h = 0.7. "
citation = "Bell et al. (2007)"
bibcode = "2007ApJ...663..834B"
name = "Star formation rates function from 24um observations."
plot_as = "line"
redshift = 0.2
h = 0.7

# Write everything
processed = ObservationalData()
processed.associate_x(
    SFR, scatter=None, comoving=True, description="Star Formation Rates"
)
processed.associate_y(Phi, scatter=None, comoving=True, description="Phi (SFR)")
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
