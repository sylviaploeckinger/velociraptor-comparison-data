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

input_filename = "../raw/Muzzin_2013.txt"
delimiter = "\t"

output_filename = "Muzzin2013.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

raw = np.loadtxt(input_filename)

# Redshift of the data point
z_lo = raw[:, 0]
z_hi = raw[:, 1]
z = (z_lo + z_hi) / 2.0

# Convert to scale factors
a = unyt.unyt_array(1 / (1 + z), "dimensionless")
a_lo = unyt.unyt_array(1 / (1 + z_hi), "dimensionless")
a_hi = unyt.unyt_array(1 / (1 + z_lo), "dimensionless")

# Read densities
rho_star = unyt.unyt_array(10 ** raw[:, 17], "Msun / Mpc**3")
rho_star_hi = unyt.unyt_array(10 ** (raw[:, 17] + raw[:, 18]), "Msun / Mpc**3")
rho_star_lo = unyt.unyt_array(10 ** (raw[:, 17] - raw[:, 19]), "Msun / Mpc**3")

# Convert from Kroupa to Chabrier IMF
rho_star = rho_star / 10 ** 0.04
rho_star_hi = rho_star / 10 ** 0.04
rho_star_low = rho_star / 10 ** 0.04

# Build scatter
x_scatter = unyt.unyt_array((a - a_lo, a_hi - a))
y_scatter = unyt.unyt_array((rho_star - rho_star_lo, rho_star_hi - rho_star))

# Meta-data
comment = (
    "Data obtained assuming a Chabrier IMF and h = 0.7. "
    "Objects with stellar mass > 10^8 Msun are considered. "
    f"h-corrected for SWIFT using cosmology: {cosmology.name}. "
)
citation = "Muzzin et al. (2013) (COSMOS)"
bibcode = "2013ApJ...777...18M"
name = "Stellar mass density obtained intagrating the GSMF from COSMOS/ultraVISTA."
redshift = z
plot_as = "points"

# Write everything
processed = ObservationalData()
processed.associate_x(a, scatter=x_scatter, comoving=True, description="Scale-factor")
processed.associate_y(
    rho_star, scatter=y_scatter, comoving=True, description="Stellar Mass Density"
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
