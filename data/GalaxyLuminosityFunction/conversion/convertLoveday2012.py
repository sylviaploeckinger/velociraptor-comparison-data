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

output_basename = "Loveday2012_"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Data taken from table 3 of the paper
bands = ["u", "g", "r", "i", "z"]
min_M = np.array([-21.0, -22.0, -23.0, -23.0, -24.0])
max_M = np.array([-10.0, -10.0, -10.0, -11.0, -12.0])
alpha = np.array([-1.21, -1.20, -1.26, -1.22, -1.18])
M_star = np.array([-18.02, -19.71, -20.73, -21.13, -21.41])
phi_star = np.array([1.96, 1.33, 0.90, 0.90, 0.90]) / 100.0

# Convert to our cosmology
min_M += 5 * np.log10(h_sim)
max_M += 5 * np.log10(h_sim)
M_star += 5 * np.log10(h_sim)
phi_star /= h_sim ** (-3)

# Meta-data
comment = f"Data h-corrected for SWIFT using cosmology: {cosmology.name}."
citation = "Loveday et al. (2012) (GAMA)"
bibcode = "2012MNRAS.420.1239L"
name = "Luminosity functions in the ugriz bands from the GAMA survey. Single-Schechter fits to the data."
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
