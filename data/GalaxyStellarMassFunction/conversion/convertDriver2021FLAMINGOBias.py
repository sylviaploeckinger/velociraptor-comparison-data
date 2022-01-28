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

output_filename = "Driver2021FLAMINGOBias.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

processed = ObservationalData()

# Table 5
log10_Mstar = np.array(
    [
        6.875,
        7.125,
        7.375,
        7.625,
        7.875,
        8.125,
        8.375,
        8.625,
        8.875,
        9.125,
        9.375,
        9.625,
        9.875,
        10.125,
        10.375,
        10.625,
        10.875,
        11.125,
        11.375,
        11.625,
    ]
)
log10_Phi = np.array(
    [
        -0.691,
        -1.084,
        -1.011,
        -1.349,
        -1.287,
        -1.544,
        -1.669,
        -1.688,
        -1.795,
        -1.886,
        -2.055,
        -2.142,
        -2.219,
        -2.274,
        -2.292,
        -2.361,
        -2.561,
        -2.922,
        -3.414,
        -4.704,
    ]
)
log10_Phi_delta = np.array(
    [
        0.176,
        0.125,
        0.071,
        0.092,
        0.079,
        0.071,
        0.045,
        0.032,
        0.024,
        0.020,
        0.014,
        0.010,
        0.009,
        0.009,
        0.009,
        0.010,
        0.013,
        0.019,
        0.032,
        0.138,
    ]
)

# Add correction to z=0 and SDSS
log10_Phi += 0.0769

#Add corrections for the FLAMINGO Biases
log10_Mstar += 0.03189831
log10_Phi += np.log10(0.99693299)

log10_Phi_plus = log10_Phi + log10_Phi_delta
log10_Phi_minus = log10_Phi - log10_Phi_delta

# Convert to SWIFT-friendly units and correct for cosmology
Mstar = (10.0 ** log10_Mstar) * unyt.Solar_Mass * (h_sim / h_obs) ** -2
Phi = (10.0 ** log10_Phi) * unyt.Mpc ** (-3) * (h_sim / h_obs) ** 3
Phi_plus = ((10.0 ** (log10_Phi_plus)) * unyt.Mpc ** (-3)) * (h_sim / h_obs) ** 3
Phi_minus = ((10.0 ** (log10_Phi_minus)) * unyt.Mpc ** (-3)) * (h_sim / h_obs) ** 3

Phi_plus = Phi_plus - Phi
Phi_minus = Phi - Phi_minus

# Meta-data
comment = (
    "Data obtained assuming a Chabrier IMF and h = 0.7. "
    f"h-corrected for SWIFT using cosmology: {cosmology.name}. "
    "Ignoring the mass bins for which GAMA is systematically incomplete."
)
citation = "Driver et al. (2021) (GAMA-DR4)"
bibcode = "private communication"
name = "GSMF from GAMA-DR4 (FLAMINGO Biases)"
redshift = 0.0
plot_as = "points"

# Write everything
processed.associate_x(
    Mstar, scatter=None, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    Phi,
    scatter=unyt.unyt_array([Phi_minus, Phi_plus]),
    comoving=True,
    description="Phi (GSMF)",
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
