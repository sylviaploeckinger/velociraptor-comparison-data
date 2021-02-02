from velociraptor.observations.objects import (
    ObservationalData,
    MultiRedshiftObservationalData,
)

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Redshifts at which to plot the data
redshifts = np.array([1.5, 2.0, 2.5])

# Valid redshift ranges for each z from above
Delta_z = 0.5 * (redshifts[1:] - redshifts[:-1])
redshifts_lower = np.append(0.10, Delta_z)
redshifts_upper = np.append(Delta_z, 0.25)

# Create the formatted version of the above array
redshift_header_info = ", ".join([f"{z:.1f}" for z in redshifts])

# Cosmology
h_sim = cosmology.h

# Meta-data
citation = "Chartab et al. (2021)"
bibcode = "2021arXiv210101706C"
name = "MOSDEF Survey: gas-phase metallicity of galaxies at 1.4 <= z <= 2.6"
plot_as = "line"
h = h_sim


name = f"Fit to the stellar mass - gas metallicity at z=[{redshift_header_info:s}]"
comment = (
    "The data is taken from Chartab+21 "
    "Median fit to galaxy stacks from MOSDEF survey. "
    "Stellar masses obtained assuming a Chabrier IMF. "
    "The metallicity is expressed as 12 + log10(O/H), in these units the solar metallicity is 8.69."
)

# Store metadata at the top level
multi_z = MultiRedshiftObservationalData()
multi_z.associate_citation(citation, bibcode)
multi_z.associate_name(name)
multi_z.associate_comment(comment)
multi_z.associate_cosmology(cosmology)
multi_z.associate_maximum_number_of_returns(1)

output_filename = "Chartab2021.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

for z, dz_lower, dz_upper in zip(redshifts, redshifts_lower, redshifts_upper):

    # Create a single observational-data instance at redshift z
    processed = ObservationalData()

    # Compute \Delta z
    redshift_lower, redshift_upper = [z - dz_lower, z + dz_upper]

    if z < 2:
        Z_asm = 6.29
        alpha = 0.21
    else:
        Z_asm = 4.84
        alpha = 0.35

    log10_M_star_min = 9.5
    log10_M_star_max = 11.5
    M_star = np.arange(log10_M_star_min, log10_M_star_max, 0.2)
    Z_gas = (Z_asm + alpha * M_star) * unyt.dimensionless  # 12 + log(O/H)
    M_star = 10 ** M_star * unyt.Solar_Mass

    processed.associate_x(
        M_star, scatter=None, comoving=True, description="Galaxy Stellar Mass"
    )
    processed.associate_y(
        Z_gas, scatter=None, comoving=True, description="Gas phase metallicity"
    )
    processed.associate_redshift(z, redshift_lower, redshift_upper)
    processed.associate_plot_as(plot_as)
    multi_z.associate_dataset(processed)

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

multi_z.write(filename=output_path)
