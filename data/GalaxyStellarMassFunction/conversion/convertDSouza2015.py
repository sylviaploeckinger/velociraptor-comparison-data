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

output_filename = "DSouza2015.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# The data is given in the form of a double Schechter fit (Table 2 of the paper)
phi_star_1 = 0.008_579 * unyt.Mpc ** (-3) * h_sim ** 3  # Mpc^-3 dex^-1
alpha_1 = -1.082
M_star_1 = (10 ** 10.615) * unyt.Solar_Mass * h_sim ** (-2)  # Msun

phi_star_2 = 0.000_355 * unyt.Mpc ** (-3) * h_sim ** 3  # Mpc^-3 dex^-1
alpha_2 = -1.120
M_star_2 = (10 ** 10.995) * unyt.Solar_Mass * h_sim ** (-2)  # Msun

# Minimal and maximal mass
M_min = (10 ** 9.5) * unyt.Solar_Mass * h_sim ** (-2)  # Msun
M_max = (10 ** 12.0) * unyt.Solar_Mass * h_sim ** (-2)  # Msun

# Create the x-data
M = np.logspace(np.log10(M_min), np.log10(M_max), 50) * unyt.Solar_Mass  # Msun

# Create the y-data (double Schechter)
Phi_Md_M = (phi_star_1 / M_star_1) * np.exp(-M / M_star_1) * (
    M / M_star_1
) ** alpha_1 + (phi_star_2 / M_star_2) * np.exp(-M / M_star_2) * (
    M / M_star_2
) ** alpha_2
Phi = Phi_Md_M * M * np.log(10)

# Meta-data
comment = (
    "Data obtained assuming a Chabrier IMF and h = 0.72. "
    f"h-corrected for SWIFT using cosmology: {cosmology.name}."
)
citation = "D'Souza et al. (2015)"
bibcode = "2015MNRAS.454.4027D"
name = "GSMF from SDSS DR7 (NYU-VAGC catalog) with flux corrections and K-corrections"
plot_as = "line"
redshift = 0.1
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(M, scatter=None, comoving=True, description="Galaxy Stellar Mass")
processed.associate_y(Phi, scatter=None, comoving=True, description="Phi (GSMF)")
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
