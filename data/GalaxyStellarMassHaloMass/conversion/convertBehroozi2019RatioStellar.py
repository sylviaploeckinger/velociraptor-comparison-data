from velociraptor.observations.objects import (
    ObservationalData,
    MultiRedshiftObservationalData,
)
from velociraptor.fitting_formulae.smhmr import behroozi_2019_raw

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Redshifts at which to plot the data
redshifts = np.array([0.0, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])

# Valid redshift ranges for each z from above
Delta_z = 0.5 * (redshifts[1:] - redshifts[:-1])
redshifts_lower = np.append(0.10, Delta_z)
redshifts_upper = np.append(Delta_z, 0.25)

# Create the formatted version of the above array
redshift_header_info = ", ".join([f"{z:.1f}" for z in redshifts])

# Halo masses (Berhoozi data were fitted in the range [10**10.5, 10**15] Msun)
M_BN98 = np.logspace(10.5, 15.0, 512)

# Cosmology
h_sim = cosmology.h

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
    "Shows the ratio between stellar mass and halo mass as a function of stellar "
    "mass. "
)
citation = "Behroozi et al. (2019)"
bibcode = "2019MNRAS.488.3143B"
name = f"Fit to the stellar mass / halo mass - stellar mass relation at z=[{redshift_header_info:s}]"
plot_as = "line"
h = h_sim

# Store metadata at the top level
multi_z = MultiRedshiftObservationalData()
multi_z.associate_citation(citation, bibcode)
multi_z.associate_name(name)
multi_z.associate_comment(comment)
multi_z.associate_cosmology(cosmology)
multi_z.associate_maximum_number_of_returns(1)

output_filename = "Behroozi2019RatioStellar.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

for z, dz_lower, dz_upper in zip(redshifts, redshifts_lower, redshifts_upper):

    # Create a single observational-data instance at redshift z
    processed = ObservationalData()

    # Stellar masses (for the given halo masses, at redshift z)
    M_star = behroozi_2019_raw(z, M_BN98)

    # Compute \Delta z
    redshift_lower, redshift_upper = [z - dz_lower, z + dz_upper]

    processed.associate_x(
        M_star * unyt.Solar_Mass,
        scatter=None,
        comoving=True,
        description="Galaxy Stellar Mass",
    )
    processed.associate_y(
        (M_star / M_BN98) * unyt.dimensionless,
        scatter=None,
        comoving=True,
        description="Galaxy Stellar Mass / Halo Mass ($M_* / M_{\\rm BN98}$)",
    )

    processed.associate_redshift(z, redshift_lower, redshift_upper)
    processed.associate_plot_as(plot_as)

    multi_z.associate_dataset(processed)

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

multi_z.write(filename=output_path)
