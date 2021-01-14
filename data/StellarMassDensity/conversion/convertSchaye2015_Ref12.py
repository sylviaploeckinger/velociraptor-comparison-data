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

input_filename = "../raw/MassEvolution_EAGLE_Ref_12.dat"
delimiter = "\t"

output_filename = "Schaye2015_Ref12.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

raw = np.loadtxt(input_filename)

# Redshift of the data points
z = raw[:, 0]

# Convert to scale factors
a = unyt.unyt_array(1 / (1 + z), "dimensionless")

# Read densities
rho_star = unyt.unyt_array(raw[:, 3], "Msun")
rho_star /= unyt.unyt_quantity(12.5, "Mpc") ** 3

# Meta-data
comment = (
    "Data obtained summing up all the star particles in the snapshots. "
    "No cosmology correction needed."
)
citation = "Schaye et al. (2015) (EAGLE-12)"
bibcode = "2015MNRAS..446..521S"
name = "Stellar mass density obtained from the Ref EAGLE model (12.5 Mpc volume)"
redshift = np.mean(z)
plot_as = "points"

# Write everything
processed = ObservationalData()
processed.associate_x(a, scatter=None, comoving=True, description="Scale-factor")
processed.associate_y(
    rho_star, scatter=None, comoving=True, description="Stellar Mass Density"
)
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(redshift)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)
processed.associate_redshift(redshift, np.min(redshift), np.max(redshift))

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)
