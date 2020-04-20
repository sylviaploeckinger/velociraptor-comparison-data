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

output_filename = "McConnell2013_Fit.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Minimal and maximal mass
M_min = (10.0 ** 8.0) * unyt.Solar_Mass  # Msun
M_max = (10.0 ** 13.0) * unyt.Solar_Mass  # Msun

# Create the x-data
M_star = np.logspace(np.log10(M_min), np.log10(M_max), 50) * unyt.Solar_Mass  # Msun

# Create the y-data (fit from the paper)
M_BH = 8.46 + 1.05 * np.log10(M_star / (1e11 * unyt.Solar_Mass))
M_BH = (10.0 ** M_BH) * unyt.Solar_Mass  # Msun

# Meta-data
comment = (
    "Data obtained assuming the total stellar mass is the same as the bulge mass. "
    "No cosmology correction needed."
)
citation = "McConnell & Ma (2013) (Fit)"
bibcode = "2013ApJ...764..184M"
name = "Fit to the black hole mass - stellar mass relation from 72 local galaxies."
plot_as = "line"
redshift = 0.0
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star, scatter=None, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(M_BH, scatter=None, comoving=True, description="Black Hole Mass")
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
