from velociraptor.observations.objects import ObservationalData
from velociraptor.fitting_formulae.smhmr import behroozi_2019_raw

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Cosmology
h_sim = cosmology.h

# Redshifts at which to plot the data
redshifts = [0.0, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
# Create the formatted version of the above array
redshift_header_info = ", ".join([f"{z:.1f}" for z in redshifts])

# Halo masses
M_BN98 = np.logspace(9, 15, 512)

# Stellar masses (for each z in redshifts)
M_star = np.array([behroozi_2019_raw(z, M_BN98) for z in redshifts])

output_filename = "Behroozi2019.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Meta-data
comment = (
    "The data is taken from https://www.peterbehroozi.com/data.html. "
    "Median fit to the raw data for centrals (i.e. excluding satellites). "
    "The stellar mass is the true stellar mass (i.e. w/o observational "
    "corrections). "
    "The halo mass is the peak halo mass that follows the Bryan & Norman (1998) "
    "spherical overdensity definition. "
    "The fitting function does not include the intrahalo light contribution to the "
    "stellar mass. "
    "Cosmology: Omega_m=0.307, Omega_lambda=0.693, h=0.678, sigma_8=0.823, "
    "n_s=0.96. "
    "Shows the stellar mass as a function halo mass."
)
citation = "Behroozi et al. (2019)"
bibcode = "2019MNRAS.488.3143B"
name = f"Fit to the stellar mass - halo mass relation at z=[{redshift_header_info:s}]"
plot_as = "line"
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_BN98 * unyt.Solar_Mass,
    scatter=None,
    comoving=True,
    description="Halo Mass ($M_{\\rm BN98}$)",
)
processed.associate_y(
    M_star * unyt.Solar_Mass,
    scatter=None,
    comoving=True,
    description="Galaxy Stellar Mass",
)
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(redshifts)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)
