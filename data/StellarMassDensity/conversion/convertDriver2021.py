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

output_filename = "Driver2021.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Value taken from page 25 of the paper
Omega_star = 2.17 * 10 ** -3 / (h_obs / h_sim)
Omega_star_p = 0.03 * 10 ** -3 / (h_obs / h_sim)
Omega_star_m = 0.03 * 10 ** -3 / (h_obs / h_sim)

# Convert to densities assuming the simulation's cosmology
rho_crit = unyt.unyt_array.from_astropy(cosmology.critical_density0).to("Msun/Mpc**3")
stellar_mass_density = unyt.unyt_array([Omega_star * rho_crit])
stellar_mass_density_p = unyt.unyt_array([Omega_star_p * rho_crit])
stellar_mass_density_m = unyt.unyt_array([Omega_star_m * rho_crit])

# Construct the error bar
y_scatter = unyt.unyt_array((stellar_mass_density_m, stellar_mass_density_p))

# Redshift of the data point
z = 0.1
a = unyt.unyt_array([1 / (1 + z)], "dimensionless")

# Meta-data
comment = (
    "Data obtained assuming a Chabrier IMF and h = 0.7. "
    f"h-corrected for SWIFT using cosmology: {cosmology.name}. "
)
citation = "Driver et al. (2021) (GAMA-DR4)"
bibcode = ""
name = "Stellar mass density obtained intagrating the GSMF from GAMA-DR4."
redshift = z
plot_as = "points"
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(a, scatter=None, comoving=True, description="Scale-factor")
processed.associate_y(
    stellar_mass_density,
    scatter=y_scatter,
    comoving=True,
    description="Stellar Mass Density",
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
