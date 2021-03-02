from velociraptor.observations.objects import (
    ObservationalData,
    MultiRedshiftObservationalData,
)

import unyt
import numpy as np
import os
import sys

# Hubble parameter used in Behroozi2019
ORIGINAL_H = 0.678


def Phi_all_galaxies():

    # Meta-data
    name = f"Fit to the galaxy stellar mass function at z=[{redshift_header_info:s}]"
    comment = (
        "The data is taken from https://www.peterbehroozi.com/data.html. "
        "The stellar mass is the observed stellar mass as defined in Behroozi et al. "
        "(2019) eq. 25. "
        "Uses the Chabrier initial mass function. "
        "GSMF is incomplete below 10**7.0 Msun at z=0 and 10**8.5 Msum at z=8. "
        "Cosmology: Omega_m=0.307, Omega_lambda=0.693, h=0.678, sigma_8=0.823, "
        "n_s=0.96. "
        "Shows the galaxy stellar mass function (number densities in comoving Mpc^-3 "
        " dex^-1 vs. stellar mass)."
    )

    # Store metadata at the top level
    multi_z = MultiRedshiftObservationalData()
    multi_z.associate_citation(citation, bibcode)
    multi_z.associate_name(name)
    multi_z.associate_comment(comment)
    multi_z.associate_cosmology(cosmology)
    multi_z.associate_maximum_number_of_returns(1)

    output_filename = "Behroozi2019_all.hdf5"
    output_directory = "../"

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    for z, dz_lower, dz_upper, a_str in zip(
        redshifts, redshifts_lower, redshifts_upper, scale_factors_str
    ):
        # Create a single observational-data instance at redshift z
        processed = ObservationalData()

        # Load raw Behroozi2019 data
        data = np.loadtxt(f"../raw/Behroozi2019_smf_a{a_str}.dat")

        # Fetch the fields we need
        log_M_star, Phi, Phi_plus, Phi_minus = (
            data[:, 0],
            data[:, 1],
            data[:, 2],
            data[:, 3],
        )

        # We don't want to plot zeros
        mask = np.where(Phi > 0.0)

        # Transform stellar mass
        M_star = (10.0 ** log_M_star) * unyt.Solar_Mass

        # Define scatter with respect to the best-fit value (16 and 84 percentiles)
        Phi_scatter = unyt.unyt_array(
            (Phi_minus[mask], Phi_plus[mask]), units=unyt.Mpc ** (-3)
        )

        # Compute \Delta z
        redshift_lower, redshift_upper = [z - dz_lower, z + dz_upper]

        processed.associate_x(
            M_star[mask],
            scatter=None,
            comoving=False,
            description="Galaxy Stellar Mass",
        )
        processed.associate_y(
            Phi[mask] * (h_sim / ORIGINAL_H) ** 3 * unyt.Mpc ** (-3),
            scatter=Phi_scatter * (h_sim / ORIGINAL_H) ** 3,
            comoving=True,
            description="Phi (GSMF)",
        )

        processed.associate_redshift(z, redshift_lower, redshift_upper)
        processed.associate_plot_as(plot_as)

        multi_z.associate_dataset(processed)

    output_path = f"{output_directory}/{output_filename}"

    if os.path.exists(output_path):
        os.remove(output_path)

    multi_z.write(filename=output_path)


def Phi_passive_galaxies():

    # Meta-data
    name = f"Fit to the quenched galaxy stellar mass function at z=[{redshift_header_info:s}]"
    comment = (
        "The data is taken from https://www.peterbehroozi.com/data.html. "
        "The stellar mass is the observed stellar mass as defined in Behroozi et al. "
        "(2019) eq. 25. "
        "The quenched fractions are defined using the standard criterion where specific"
        " star formation rate < 1e-11 yr^-1. "
        "Uses the Chabrier initial mass function. "
        "GSMF is incomplete below 10**7.0 Msun at z=0 and 10**8.5 Msum at z=8. "
        "Cosmology: Omega_m=0.307, Omega_lambda=0.693, h=0.678, sigma_8=0.823, "
        "n_s=0.96. "
        "Shows the quenched galaxy stellar mass function (number densities in comoving"
        " Mpc^-3 dex^-1 vs. stellar mass)."
    )

    # Store metadata at the top level
    multi_z = MultiRedshiftObservationalData()
    multi_z.associate_citation(citation, bibcode)
    multi_z.associate_name(name)
    multi_z.associate_comment(comment)
    multi_z.associate_cosmology(cosmology)
    multi_z.associate_maximum_number_of_returns(1)

    output_filename = "Behroozi2019_passive.hdf5"
    output_directory = "../"

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    for z, dz_lower, dz_upper, a_str in zip(
        redshifts, redshifts_lower, redshifts_upper, scale_factors_str
    ):
        # Create a single observational-data instance at redshift z
        processed = ObservationalData()

        # Load raw Behroozi2019 data
        data = np.loadtxt(f"../raw/Behroozi2019_smf_a{a_str}.dat")

        # Fetch the fields we need
        log_M_star, Phi, Phi_plus, Phi_minus = (
            data[:, 0],
            data[:, 7],
            data[:, 8],
            data[:, 9],
        )

        # We don't want to plot zeros
        mask = np.where(Phi > 0.0)

        # Transform stellar mass
        M_star = (10.0 ** log_M_star) * unyt.Solar_Mass

        # Define scatter with respect to the best-fit value (16 and 84 percentiles)
        Phi_scatter = unyt.unyt_array(
            (Phi_minus[mask], Phi_plus[mask]), units=unyt.Mpc ** (-3)
        )

        # Compute \Delta z
        redshift_lower, redshift_upper = [z - dz_lower, z + dz_upper]

        processed.associate_x(
            M_star[mask],
            scatter=None,
            comoving=False,
            description="Galaxy Stellar Mass",
        )
        processed.associate_y(
            Phi[mask] * (h_sim / ORIGINAL_H) ** 3 * unyt.Mpc ** (-3),
            scatter=Phi_scatter * (h_sim / ORIGINAL_H) ** 3,
            comoving=True,
            description="Phi (GSMF)",
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

# Cosmology
h_sim = cosmology.h

# Scale factors at which to plot the data
scale_factors_str = [
    "1.002310",
    "0.835248",
    "0.658060",
    "0.501122",
    "0.399872",
    "0.328997",
    "0.283435",
    "0.247998",
    "0.222685",
    "0.202435",
    "0.182185",
    "0.166998",
]
scale_factors = np.array([float(a) for a in scale_factors_str])

redshifts = 1.0 / scale_factors - 1.0

# Valid redshift ranges for each z from above
Delta_z = 0.5 * (redshifts[1:] - redshifts[:-1])
redshifts_lower = np.append(0.10, Delta_z)
redshifts_upper = np.append(Delta_z, 0.25)

# Create the formatted version of the above array
redshift_header_info = ", ".join([f"{z:.1f}" for z in redshifts])

# Metadata
citation = "Behroozi et al. (2019)"
bibcode = "2019MNRAS.488.3143B"
plot_as = "line"

# Generate, format and save the Behroozi2019 data
Phi_all_galaxies()
Phi_passive_galaxies()
