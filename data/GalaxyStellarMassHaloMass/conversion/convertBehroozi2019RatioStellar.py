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
redshifts = [0.0, 0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0]

# Create and save data for each z in redshifts
for z in redshifts:

    output_filename = f"Behroozi2019RatioStellar_{f'z{z:07.3f}'.replace('.', 'p')}.hdf5"
    output_directory = "../"

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    M_BN98 = np.logspace(9, 15, 512)

    # Note that the function below actually takes M_{BN98, peak}, and we are
    # ignoring this fact
    M_star = behroozi_2019_raw(z, M_BN98)

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
    name = "Fit to the stellar mass - halo mass relation at z={:.1f}".format(z)
    plot_as = "line"
    h = h_sim

    # Write everything
    processed = ObservationalData()
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
    processed.associate_citation(citation, bibcode)
    processed.associate_name(name)
    processed.associate_comment(comment)
    processed.associate_redshift(z)
    processed.associate_plot_as(plot_as)
    processed.associate_cosmology(cosmology)

    output_path = f"{output_directory}/{output_filename}"

    if os.path.exists(output_path):
        os.remove(output_path)

    processed.write(filename=output_path)
