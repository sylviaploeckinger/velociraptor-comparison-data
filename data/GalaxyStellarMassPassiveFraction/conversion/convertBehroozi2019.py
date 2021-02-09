from velociraptor.observations.objects import (
    ObservationalData,
    MultiRedshiftObservationalData,
)

import unyt
import numpy as np
import os
import sys


def passive_fractions():

    # Meta-data
    name = f"Fit to the passive fraction - stellar mass at z=[{redshift_header_info:s}]"
    comment = (
        "The data is taken from https://www.peterbehroozi.com/data.html. "
        "The quenched fractions are defined using the standard criterion where "
        "specific star formation rate < 1e-11 yr^-1. "
        "The stellar mass is the observed stellar mass as defined in Behroozi et al. "
        "(2019) eq. 25. "
        "Uses the Chabrier initial mass function. "
        "The passive fractions are provided for the best-fitting model. In certain "
        "rare cases, the best model may lie outside the 16th-84th percentile range. "
        "Cosmology: Omega_m=0.307, Omega_lambda=0.693, h=0.678, sigma_8=0.823, "
        "n_s=0.96. "
        "Shows the passive fraction of galaxies versus galaxy stellar mass."
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

    for z, dz_lower, dz_upper, a_str in zip(
        redshifts, redshifts_lower, redshifts_upper, scale_factors_str
    ):
        # Create a single observational-data instance at redshift z
        processed = ObservationalData()

        # Load raw Behroozi2019 data
        data = np.loadtxt(f"../raw/Behroozi2019_qf_a{a_str}.dat")

        # Fetch the fields we need
        log_M_star, QF, QF_plus, QF_minus = (
            data[:, 0],
            data[:, 4],
            data[:, 5],
            data[:, 6],
        )

        # We don't want to plot zeros
        mask = np.where(QF > 0.0)

        # Transform stellar mass
        M_star = (10.0 ** log_M_star) * unyt.Solar_Mass

        # Define scatter with respect to the best-fit value (16 and 84 percentiles)
        QF_scatter = unyt.unyt_array(
            (QF_minus[mask], QF_plus[mask]), units="dimensionless"
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
            QF[mask] * unyt.dimensionless,
            scatter=QF_scatter,
            comoving=False,
            description="Passive Fraction",
        )

        processed.associate_redshift(z, redshift_lower, redshift_upper)
        processed.associate_plot_as(plot_as)

        multi_z.associate_dataset(processed)

    output_path = f"{output_directory}/{output_filename}"

    if os.path.exists(output_path):
        os.remove(output_path)

    multi_z.write(filename=output_path)


def passive_fractions_centrals():

    # Meta-data
    name = (
        "Fit to the passive fraction - stellar mass (centrals) "
        f"at z=[{redshift_header_info:s}]"
    )
    comment = (
        "The data is taken from https://www.peterbehroozi.com/data.html. "
        "The quenched fractions are defined using the standard criterion where "
        "specific star formation rate < 1e-11 yr^-1. "
        "The stellar mass is the observed stellar mass as defined in Behroozi et al. "
        "(2019) eq. 25. "
        "Uses the Chabrier initial mass function. "
        "The passive fractions are given by the 50th percentile of the posterior "
        "distribution of the fitting model. "
        "Cosmology: Omega_m=0.307, Omega_lambda=0.693, h=0.678, sigma_8=0.823, "
        "n_s=0.96. "
        "Shows the passive fraction of centrals versus galaxy stellar mass."
    )

    # Store metadata at the top level
    multi_z = MultiRedshiftObservationalData()
    multi_z.associate_citation(citation, bibcode)
    multi_z.associate_name(name)
    multi_z.associate_comment(comment)
    multi_z.associate_cosmology(cosmology)
    multi_z.associate_maximum_number_of_returns(1)

    output_filename = "Behroozi2019_centrals.hdf5"
    output_directory = "../"

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    for z, dz_lower, dz_upper, a_str in zip(
        redshifts, redshifts_lower, redshifts_upper, scale_factors_str
    ):
        # Create a single observational-data instance at redshift z
        processed = ObservationalData()

        # Load raw Behroozi2019 data
        data = np.loadtxt(f"../raw/Behroozi2019_qf_groupstats_a{a_str}.dat")

        # Fetch the fields we need
        log_M_star, QF, QF_plus, QF_minus = (
            data[:, 0],
            data[:, 1],
            data[:, 2],
            data[:, 3],
        )

        # We don't want to plot zeros
        mask = np.where(QF > 0.0)

        # Transform stellar mass
        M_star = (10.0 ** log_M_star) * unyt.Solar_Mass

        # Define scatter with respect to the best-fit value (16 and 84 percentiles)
        QF_scatter = unyt.unyt_array(
            (QF_minus[mask], QF_plus[mask]), units="dimensionless"
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
            QF[mask] * unyt.dimensionless,
            scatter=QF_scatter,
            comoving=False,
            description="Passive Fraction (centrals)",
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
passive_fractions()
passive_fractions_centrals()
