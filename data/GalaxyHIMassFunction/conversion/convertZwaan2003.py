from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

output_filename = "Zwaan2003.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# The data is given in the form of a Schechter fit (Fig. 4 of the paper)
phi_star = unyt.unyt_quantity(0.0086, units="Mpc**-3")
alpha = -1.30
M_star = unyt.unyt_quantity(10 ** 9.79, units="Msun")

# Minimal and maximal mass
M_min = unyt.unyt_quantity(10 ** 6.7, units="Msun")
M_max = unyt.unyt_quantity(10 ** 10.7, units="Msun")

# Create the x-data
M_HI = unyt.unyt_array(np.logspace(np.log10(M_min), np.log10(M_max), 100), units="Msun")

# Create the y-data (Schechter function)
Phi_HI = (phi_star / M_star) * np.exp(-M_HI / M_star) * (M_HI / M_star) ** alpha
Phi_HI = Phi_HI * M_HI * np.log(10)

# Meta-data
comment = "Best-fit Schechter function to the data from 1000 galaxies"
citation = "Zwaan et al. (2003) (HIPASS)"
bibcode = "2003AJ....125.2842Z"
name = "HI mass function from the HIPASS survey"
plot_as = "line"
redshift = 0.0
h = 0.7

# Write everything
processed = ObservationalData()
processed.associate_x(M_HI, scatter=None, comoving=True, description="HI masses")
processed.associate_y(Phi_HI, scatter=None, comoving=True, description="Phi (M_HI)")
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
