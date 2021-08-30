from velociraptor.observations.objects import (
    ObservationalData,
    MultiRedshiftObservationalData,
)

import unyt
import numpy as np
import os
import sys


def cosmic_star_formation_history_observed():

    # Meta-data
    name = f"Fit to observed star formation rate history from Behroozi et al. (2019)"
    comment = (
        "The data is taken from https://www.peterbehroozi.com/data.html. "
        "The cosmic star formation rate is the observed one. This means that "
        "systematic uncertainties in stellar masses and SFRs that arise from "
        "modeling assumptions for stellar population synthesis, dust, "
        "metallicity, and star formation history are accounted for using a "
        "redshift-dependent correction factor [See eqs. 25 and 26 in Berhoozi "
        "et al. (2019)]. "
        "Uses the Chabrier initial mass function. "
        "The scatter shows the 16th-84th percentile range from the posterior. "
        "Cosmology: Omega_m=0.307, Omega_lambda=0.693, h=0.678, sigma_8=0.823, "
        "n_s=0.96. "
        "Shows total observed cosmic star formation rate (Msun/yr/Mpc^3) for "
        "the best-fitting model from Behroozi et al. (2019)"
    )

    citation = "Behroozi et al. (2019) [Observed]"
    bibcode = "2019MNRAS.488.3143B"
    plot_as = "line"
    output_filename = "Behroozi2019_observed.hdf5"
    output_directory = "../"

    # Create observational data instance
    processed = ObservationalData()
    processed.associate_citation(citation, bibcode)
    processed.associate_name(name)
    processed.associate_comment(comment)
    processed.associate_cosmology(cosmology)

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    # Load raw Behroozi2019 data
    data = np.loadtxt(f"../raw/Behroozi2019_csfrs.dat")

    # Fetch the fields we need
    scale_factor = unyt.unyt_array(data[:, 0], units="dimensionless")
    SFR, SFR_plus, SFR_minus = data[:, 1], data[:, 3], data[:, 2]

    # Define scatter with respect to the best-fit value (16 and 84 percentiles)
    SFR_scatter = unyt.unyt_array((SFR_minus, SFR_plus), units="Msun/yr/Mpc**3")

    redshift, redshift_lower, redshift_upper = 5.0, 0.0, 10.0

    processed.associate_x(
        scale_factor,
        scatter=None,
        comoving=False,
        description="Cosmic scale factor",
    )
    processed.associate_y(
        SFR * SFR_scatter.units,
        scatter=SFR_scatter,
        comoving=False,
        description="Cosmic average star formation rate density",
    )

    processed.associate_redshift(redshift, redshift_lower, redshift_upper)
    processed.associate_plot_as(plot_as)

    output_path = f"{output_directory}/{output_filename}"

    if os.path.exists(output_path):
        os.remove(output_path)

    processed.write(filename=output_path)


def cosmic_star_formation_history_true():

    # Meta-data
    name = f"Fit to true star formation rate history from Behroozi et al. (2019)"
    comment = (
        "The data is taken from https://www.peterbehroozi.com/data.html. "
        "The cosmic star formation rate is the true one. This means that no"
        "observational systematics (e.g., due to dust) are accounted for. "
        "Uses the Chabrier initial mass function. "
        "The scatter shows the 16th-84th percentile range from the posterior. "
        "Cosmology: Omega_m=0.307, Omega_lambda=0.693, h=0.678, sigma_8=0.823, "
        "n_s=0.96. "
        "Shows total true cosmic star formation rate (Msun/yr/Mpc^3) for "
        "the best-fitting model from Behroozi et al. (2019)"
    )

    citation = "Behroozi et al. (2019) [True]"
    bibcode = "2019MNRAS.488.3143B"
    plot_as = "line"
    output_filename = "Behroozi2019_true.hdf5"
    output_directory = "../"

    # Create observational data instance
    processed = ObservationalData()
    processed.associate_citation(citation, bibcode)
    processed.associate_name(name)
    processed.associate_comment(comment)
    processed.associate_cosmology(cosmology)

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    # Load raw Behroozi2019 data
    data = np.loadtxt(f"../raw/Behroozi2019_csfrs.dat")

    # Fetch the fields we need
    scale_factor = unyt.unyt_array(data[:, 0], units="dimensionless")
    SFR, SFR_plus, SFR_minus = data[:, 7], data[:, 9], data[:, 8]

    # Define scatter with respect to the best-fit value (16 and 84 percentiles)
    SFR_scatter = unyt.unyt_array((SFR_minus, SFR_plus), units="Msun/yr/Mpc**3")

    redshift, redshift_lower, redshift_upper = 5.0, 0.0, 10.0

    processed.associate_x(
        scale_factor,
        scatter=None,
        comoving=False,
        description="Cosmic scale factor",
    )
    processed.associate_y(
        SFR * SFR_scatter.units,
        scatter=SFR_scatter,
        comoving=False,
        description="Cosmic average star formation rate density",
    )

    processed.associate_redshift(redshift, redshift_lower, redshift_upper)
    processed.associate_plot_as(plot_as)

    output_path = f"{output_directory}/{output_filename}"

    if os.path.exists(output_path):
        os.remove(output_path)

    processed.write(filename=output_path)

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Generate, format and save the Behroozi2019 data
cosmic_star_formation_history_observed()
cosmic_star_formation_history_true()
