from velociraptor.observations.objects import (
    ObservationalData,
    MultiRedshiftObservationalData,
)

import unyt
import numpy as np
import os
import sys
from typing import Tuple, Any


def moster_2018_ratios(
    z: float, Mhalo: Any, path_to_param_file: str
) -> Tuple[Any, Any, Any]:
    """
    Stellar-to-halo mass ratio vs. halo mass from Moster +2018.

    The data is taken from MNRAS, Volume 477, Issue 2, p.1822-1852,
    Moster, Benjamin P. ; Naab, Thorsten ; White, Simon D. M, table 8, all centrals

    This function is a median fit to the raw data for centrals (i.e. excluding
    satellites)

    The halo mass is the peak halo mass that follows the Bryan & Norman (1998)
    spherical overdensity definition

    The function returns the best-fit stellar-to-halo mass ratio as well as
    its 68% confidence interval range

    The function must be provided with a .txt file containing the fitting parameters
    taken directly from table 8 in Moster +2018

    Provided by Evgenii Chaikin.

    Parameters
    ----------
    z: float
        Redshift at which SHMH relation is computed
    Mhalo: Any
        Halo mass
    path_to_param_file: str
        Path to the file containing the fitting parameters

    Returns
    -------
    output: Tuple[Any, Any, Any]
        A three-length tuple containing the best-fit stellar-to-halo mass ratio,
        the 84-percentile of the stellar-to-halo mass ratio,
        the 16-percentile of the stellar-to-halo mass ratio
    """

    # Universal baryon fraction (Omega_b / Omega_m) from Moster +2018
    f_baryon = 0.156

    # Load the fitting parameters
    params_list = np.loadtxt(path_to_param_file)

    # Expected names of the parameters from the .txt file
    param_names = "z M1 epsilon_N beta gamma M_sigma sigma_0 alpha".split(" ")

    # Convert everything into a dict for better readability
    fitting_params = {name: params_list[:, i] for i, name in enumerate(param_names)}

    # Find z in the loaded data that is closest to the provided z
    z_idx = np.argmin(np.abs(z - fitting_params["z"]))

    # Sanity check
    assert np.abs(z - fitting_params["z"][z_idx]) < 1e-5, (
        f"Provided redshift (z={z}) cannot be found "
        f"in the data. The closest one is "
        f'z = {fitting_params["z"][z_idx]}'
    )

    # Fit to the integrated efficiency (Eq. 5 in Moster +2018)
    epsilon = (
        2.0
        * fitting_params["epsilon_N"][z_idx]
        * np.power(
            np.power(
                Mhalo / 10.0 ** fitting_params["M1"][z_idx],
                -fitting_params["beta"][z_idx],
            )
            + np.power(
                Mhalo / 10.0 ** fitting_params["M1"][z_idx],
                fitting_params["gamma"][z_idx],
            ),
            -1.0,
        )
    )

    # Fit to the logarithmic scatter (in dex) (Eq. 25 in Moster +2018)
    sigma = fitting_params["sigma_0"][z_idx] + np.log10(
        np.power(
            Mhalo / 10.0 ** fitting_params["M_sigma"][z_idx],
            -fitting_params["alpha"][z_idx],
        )
        + 1.0
    )

    # Compute stellar-to-halo mass ratios
    MstarMhalo = unyt.unyt_array(epsilon * f_baryon, units="dimensionless")

    # Compute one-sigma deviations from the median
    MstarMhalo_upper = unyt.unyt_array(
        epsilon * np.power(10, sigma) * f_baryon, units="dimensionless"
    )
    MstarMhalo_lower = unyt.unyt_array(
        epsilon * np.power(10, -sigma) * f_baryon, units="dimensionless"
    )

    return MstarMhalo, MstarMhalo_upper, MstarMhalo_lower


# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Redshifts at which to plot the data (same as those in Fig. 12 in Moster +2018)
redshifts = np.array([0.1, 0.5, 1.0, 2.0, 4.0, 8.0])

# Valid redshift ranges for each z from above
Delta_z = 0.5 * (redshifts[1:] - redshifts[:-1])
redshifts_lower = np.append(0.10, Delta_z)
redshifts_upper = np.append(Delta_z, 2.0)

# Create the formatted version of the above array
redshift_header_info = ", ".join([f"{z:.1f}" for z in redshifts])

# The range of halo masses used in Moster +2018 (Msun)
M_BN98 = np.logspace(10.5, 15.0, 35)

# Cosmology
h_sim = cosmology.h

for x_axis_variable in ["halo mass", "stellar mass"]:

    # Meta-data
    name = f"Fit to the stellar mass / halo mass - {x_axis_variable:s} relation at z=[{redshift_header_info:s}]"
    comment = (
        "The data is taken from MNRAS, Volume 477, Issue 2, p.1822-1852, "
        "Moster, Benjamin P. ; Naab, Thorsten ; White, Simon D. M, table 8, all centrals "
        "Median fit to the data for centrals (i.e. excluding satellites). "
        "The halo mass is the peak halo mass that follows the Bryan & Norman (1998) "
        "spherical overdensity definition. "
        "The fitting function does not include intra-cluster mass contribution to the "
        "stellar mass. "
        "Cosmology: Omega_m=0.308, Omega_lambda=0.692, h=0.6781, sigma_8=0.8149, "
        "n_s=0.9677, Omega_b=0.0484. "
        f"Shows the ratio between stellar mass and halo mass as a function of {x_axis_variable:s}."
    )

    citation = "Moster et al. (2018)"
    bibcode = "2018MNRAS.477.1822M"
    plot_as = "line"
    h = h_sim

    # Store metadata at the top level
    multi_z = MultiRedshiftObservationalData()
    multi_z.associate_citation(citation, bibcode)
    multi_z.associate_name(name)
    multi_z.associate_comment(comment)
    multi_z.associate_cosmology(cosmology)
    multi_z.associate_maximum_number_of_returns(1)

    if x_axis_variable == "halo mass":
        output_filename = "Moster2018Ratio.hdf5"
    elif x_axis_variable == "stellar mass":
        output_filename = "Moster2018RatioStellar.hdf5"
    else:
        raise ValueError(f"x_axis_variable has incorrect value ({x_axis_variable:s})")

    output_directory = "../"

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    for z, dz_lower, dz_upper in zip(redshifts, redshifts_lower, redshifts_upper):
        # Create a single observational-data instance at redshift z
        processed = ObservationalData()

        # Stellar-to-halo mass ratios (for the given halo masses, at redshift z)
        MstarMhalo, MstarMhalo_84, MstarMhalo_16 = moster_2018_ratios(
            z, M_BN98, "../raw/Moster2018.txt"
        )

        # Compute \Delta z
        redshift_lower, redshift_upper = [z - dz_lower, z + dz_upper]

        # Define scatter
        y_scatter = unyt.unyt_array(
            (MstarMhalo - MstarMhalo_16, MstarMhalo_84 - MstarMhalo)
        )

        if x_axis_variable == "halo mass":
            processed.associate_x(
                M_BN98 * unyt.Solar_Mass,
                scatter=None,
                comoving=True,
                description="Halo Mass ($M_{\\rm BN98}$)",
            )
        elif x_axis_variable == "stellar mass":
            M_star = MstarMhalo * (M_BN98 * unyt.Solar_Mass)
            processed.associate_x(
                M_star,
                scatter=None,
                comoving=True,
                description="Galaxy Stellar Mass",
            )
        else:
            raise ValueError(
                f"x_axis_variable has incorrect value ({x_axis_variable:s})"
            )

        processed.associate_y(
            MstarMhalo,
            scatter=y_scatter,
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
