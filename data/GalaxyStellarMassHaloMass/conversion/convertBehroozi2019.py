from velociraptor.observations.objects import (
    ObservationalData,
    MultiRedshiftObservationalData,
)

import unyt
import numpy as np
import os
import sys
from typing import Tuple, Any


def behroozi_2019_raw_with_uncertainties(
    z: float, Mhalo: Any, path_to_param_file: str
) -> Tuple[Any, Any, Any]:
    """
    Stellar mass-halo mass relation from Behroozi +2019.

    The data is taken from https://www.peterbehroozi.com/data.html

    This function is a median fit to the raw data for centrals (i.e. excluding
    satellites)

    The stellar mass is the true stellar mass (i.e. w/o observational corrections)

    The fitting function does not include the intrahalo light contribution to the
    stellar mass

    The halo mass is the peak halo mass that follows the Bryan & Norman (1998)
    spherical overdensity definition

    The function the best-fit smellar-to-halo mass ratios as well as their
    68% confidence interval range

    The function must be provided with a .txt file with the fitting parameters taken
    directly from https://www.peterbehroozi.com/data.html

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
        A three-length tuple containing the best-fit stellar mass,
        the 84-percentile of the stellar mass, the 16-percentile of the stellar mass
    """

    Mhalo_log = np.log10(Mhalo)

    # Load the fitting parameters:
    params_list = np.loadtxt(path_to_param_file)

    # Ensure the file has correct dimensions
    assert np.shape(params_list)[1] >= 20

    param_names = (
        "EFF_0 EFF_0_A EFF_0_A2 EFF_0_Z M_1 M_1_A "
        "M_1_A2 M_1_Z ALPHA ALPHA_A ALPHA_A2 ALPHA_Z "
        "BETA BETA_A BETA_Z DELTA GAMMA GAMMA_A GAMMA_Z CHI2".split(" ")
    )

    a = 1.0 / (1.0 + z)
    a1 = a - 1.0
    lna = np.log(a)

    try:
        log_mstar_all = np.zeros((len(params_list), len(Mhalo)))
    # In case Mhalo does not have __len__ method
    except TypeError:
        log_mstar_all = np.zeros((len(params_list), 1))

    for count, params in enumerate(params_list):

        params = dict(zip(param_names, params))

        zparams = {}

        zparams["m_1"] = (
            params["M_1"]
            + a1 * params["M_1_A"]
            - lna * params["M_1_A2"]
            + z * params["M_1_Z"]
        )
        zparams["sm_0"] = (
            zparams["m_1"]
            + params["EFF_0"]
            + a1 * params["EFF_0_A"]
            - lna * params["EFF_0_A2"]
            + z * params["EFF_0_Z"]
        )
        zparams["alpha"] = (
            params["ALPHA"]
            + a1 * params["ALPHA_A"]
            - lna * params["ALPHA_A2"]
            + z * params["ALPHA_Z"]
        )
        zparams["beta"] = params["BETA"] + a1 * params["BETA_A"] + z * params["BETA_Z"]
        zparams["delta"] = params["DELTA"]
        zparams["gamma"] = 10 ** (
            params["GAMMA"] + a1 * params["GAMMA_A"] + z * params["GAMMA_Z"]
        )

        dm = Mhalo_log - zparams["m_1"]
        dm2 = dm / zparams["delta"]
        logmstar = (
            zparams["sm_0"]
            - np.log10(10 ** (-zparams["alpha"] * dm) + 10 ** (-zparams["beta"] * dm))
            + zparams["gamma"] * np.exp(-0.5 * (dm2 * dm2))
        )

        log_mstar_all[count, :] = logmstar

    # The best-fit stellar mass
    log_mstar_best = np.copy(log_mstar_all[0, :])

    for i in range(np.shape(log_mstar_all)[1]):
        log_mstar_all[:, i] = sorted(log_mstar_all[:, i], key=lambda x: x)

    log_mstar_84 = log_mstar_all[int((1 + 0.6827) * len(log_mstar_all) / 2.0)]
    log_mstar_16 = log_mstar_all[int((1 - 0.6827) * len(log_mstar_all) / 2.0)]

    return 10.0 ** log_mstar_best, 10.0 ** log_mstar_84, 10.0 ** log_mstar_16


def StellarMass_vs_HaloMass():

    name = (
        f"Fit to the stellar mass - halo mass relation at z=[{redshift_header_info:s}]"
    )
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

    # Store metadata at the top level
    multi_z = MultiRedshiftObservationalData()
    multi_z.associate_citation(citation, bibcode)
    multi_z.associate_name(name)
    multi_z.associate_comment(comment)
    multi_z.associate_cosmology(cosmology)
    multi_z.associate_maximum_number_of_returns(1)

    output_filename = "Behroozi2019.hdf5"
    output_directory = "../"

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    for z, dz_lower, dz_upper in zip(redshifts, redshifts_lower, redshifts_upper):

        # Create a single observational-data instance at redshift z
        processed = ObservationalData()

        # Stellar masses (for the given halo masses, at redshift z)
        M_star, M_84, M_16 = behroozi_2019_raw_with_uncertainties(
            z, M_BN98, "../raw/Behroozi2019_params_smhm_true_med_cen_params.txt"
        )

        # Compute \Delta z
        redshift_lower, redshift_upper = [z - dz_lower, z + dz_upper]

        # Define scatter
        y_scatter = unyt.unyt_array((M_star - M_16, M_84 - M_star))

        processed.associate_x(
            M_BN98 * unyt.Solar_Mass,
            scatter=None,
            comoving=True,
            description="Halo Mass ($M_{\\rm BN98}$)",
        )
        processed.associate_y(
            M_star * unyt.Solar_Mass,
            scatter=y_scatter * unyt.Solar_Mass,
            comoving=True,
            description="Galaxy Stellar Mass",
        )

        processed.associate_redshift(z, redshift_lower, redshift_upper)
        processed.associate_plot_as(plot_as)

        multi_z.associate_dataset(processed)

    output_path = f"{output_directory}/{output_filename}"

    if os.path.exists(output_path):
        os.remove(output_path)

    multi_z.write(filename=output_path)


def StellarMassHaloMassRatios_vs_HaloMass():

    name = f"Fit to the stellar mass / halo mass - halo mass relation at z=[{redshift_header_info:s}]"
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
        "Shows the ratio between stellar mass and halo mass as a function of halo mass."
    )

    # Store metadata at the top level
    multi_z = MultiRedshiftObservationalData()
    multi_z.associate_citation(citation, bibcode)
    multi_z.associate_name(name)
    multi_z.associate_comment(comment)
    multi_z.associate_cosmology(cosmology)
    multi_z.associate_maximum_number_of_returns(1)

    output_filename = "Behroozi2019Ratio.hdf5"
    output_directory = "../"

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    for z, dz_lower, dz_upper in zip(redshifts, redshifts_lower, redshifts_upper):
        # Create a single observational-data instance at redshift z
        processed = ObservationalData()

        # Stellar masses (for the given halo masses, at redshift z)
        # Stellar masses (for the given halo masses, at redshift z)
        M_star, M_84, M_16 = behroozi_2019_raw_with_uncertainties(
            z, M_BN98, "../raw/Behroozi2019_params_smhm_true_med_cen_params.txt"
        )

        # Compute \Delta z
        redshift_lower, redshift_upper = [z - dz_lower, z + dz_upper]

        # Define scatter
        y_scatter = unyt.unyt_array(
            ((M_star - M_16) / M_BN98, (M_84 - M_star) / M_BN98)
        )

        processed.associate_x(
            M_BN98 * unyt.Solar_Mass,
            scatter=None,
            comoving=True,
            description="Halo Mass ($M_{\\rm BN98}$)",
        )
        processed.associate_y(
            (M_star / M_BN98) * unyt.dimensionless,
            scatter=y_scatter * unyt.dimensionless,
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


def StellarMassHaloMassRatios_vs_StellarMass():

    name = f"Fit to the stellar mass / halo mass - stellar mass relation at z=[{redshift_header_info:s}]"
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
        # Stellar masses (for the given halo masses, at redshift z)
        M_star, M_84, M_16 = behroozi_2019_raw_with_uncertainties(
            z, M_BN98, "../raw/Behroozi2019_params_smhm_true_med_cen_params.txt"
        )

        # Compute \Delta z
        redshift_lower, redshift_upper = [z - dz_lower, z + dz_upper]

        # Define scatter
        y_scatter = unyt.unyt_array(
            ((M_star - M_16) / M_BN98, (M_84 - M_star) / M_BN98)
        )

        processed.associate_x(
            M_star * unyt.Solar_Mass,
            scatter=None,
            comoving=True,
            description="Galaxy Stellar Mass",
        )
        processed.associate_y(
            (M_star / M_BN98) * unyt.dimensionless,
            scatter=y_scatter * unyt.dimensionless,
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
M_BN98 = np.logspace(10.5, 15.0, 30)

# Cosmology
h_sim = cosmology.h

# Meta-data
citation = "Behroozi et al. (2019)"
bibcode = "2019MNRAS.488.3143B"
plot_as = "line"
h = h_sim

# Create, format, and save the data for three different cases
StellarMass_vs_HaloMass()
StellarMassHaloMassRatios_vs_HaloMass()
StellarMassHaloMassRatios_vs_StellarMass()
