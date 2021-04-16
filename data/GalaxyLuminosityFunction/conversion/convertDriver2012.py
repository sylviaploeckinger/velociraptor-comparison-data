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

output_basename = "Driver2012_"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Data taken from table 4 of the paper
bands = ["u", "g", "r", "i", "z", "Y", "J", "H", "K"]
min_M = np.ones(len(bands)) * (-24.0)
max_M = np.ones(len(bands)) * (-12.0)
alpha = np.array([-1.03, -1.10, -1.12, -1.17, -1.14, -1.12, -1.10, -1.07, -1.03])
M_star = np.array(
    [-18.60, -20.09, -20.86, -21.30, -21.52, -21.63, -21.74, -21.99, -21.63]
)
phi_star = np.array([2.03, 1.47, 1.24, 1.00, 1.02, 0.98, 0.97, 1.03, 1.10]) / 100.0

# Convert to our cosmology
min_M += 5 * np.log10(h_sim)
max_M += 5 * np.log10(h_sim)
M_star += 5 * np.log10(h_sim)
phi_star /= h_sim ** (-3)

# Meta-data
comment = f"Data h-corrected for SWIFT using cosmology: {cosmology.name}."
citation = "Driver et al. (2012) (GAMA)"
bibcode = "2012MNRAS.427.3244D"
name = "Luminosity functions in the ugrizYJHK bands from the GAMA survey. Single-Schechter fits to the data."
plot_as = "line"
redshift = 0.1
h = h_sim

for i in range(len(bands)):

    M = np.linspace(min_M[i], max_M[i], 256)
    arg = 0.4 * (M_star[i] - M)
    Phi = (
        unyt.Mpc ** (-3)
        * 0.4
        * np.log(10)
        * phi_star[i]
        * ((10.0 ** arg) ** (1.0 + alpha[i]))
        * np.exp(-10.0 ** arg)
    )
    M = unyt.unyt_array(M, "dimensionless")

    # Write everything
    processed = ObservationalData()
    processed.associate_x(
        M, scatter=None, comoving=True, description="Magnitudes %s-band" % bands[i]
    )
    processed.associate_y(Phi, scatter=None, comoving=True, description="Phi (M)")
    processed.associate_citation(citation, bibcode)
    processed.associate_name(name)
    processed.associate_comment(comment)
    processed.associate_redshift(redshift)
    processed.associate_plot_as(plot_as)
    processed.associate_cosmology(cosmology)

    output_filename = output_basename + bands[i] + ".hdf5"
    output_path = f"{output_directory}/{output_filename}"

    if os.path.exists(output_path):
        os.remove(output_path)

    processed.write(filename=output_path)
