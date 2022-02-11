from velociraptor.observations.objects import (
    ObservationalData,
    MultiRedshiftObservationalData,
)

import unyt
import numpy as np
import os
import sys


def cosmic_star_formation_history_novak():

    # Meta-data
    name = f"Star formation rate density from Novak et al. (2017)"
    comment = (
        "based on JVLA COSMOS radio observations at 3 GHz "
        "cosmology is H0=70, OmegaM=0.3, OmegaL=0.7 "
        "Uses a Chabrier IMF"
    )

    citation = "Novak et al. (2017)"
    bibcode = "2017A&A...602A...5N"
    plot_as = "points"
    output_filename = "Novak2017.hdf5"
    output_directory = "../"

    # Create observational data instance
    processed = ObservationalData()
    processed.associate_citation(citation, bibcode)
    processed.associate_name(name)
    processed.associate_comment(comment)
    processed.associate_cosmology(cosmology)

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    # Load raw Novak2017 data
    data = np.loadtxt(f"../raw/sfr_novak2017.dat")

    # Fetch the fields we need
    z, delta_z_minus, delta_z_plus = data[:, 0], data[:, 1], data[:, 2]
    SFR, delta_SFR_minus, delta_SFR_plus = data[:, 3], data[:, 4], data[:, 5]

    a = 1.0 / (1.0 + z)
    # division turns large into small, so minus <-> plus
    a_minus = 1.0 / (1.0 + z + delta_z_plus)
    a_plus = 1.0 / (1.0 + z + delta_z_minus)
    delta_a_minus = a - a_minus
    delta_a_plus = a_plus - a

    a_bin = unyt.unyt_array(a, units="dimensionless")
    a_scatter = unyt.unyt_array((delta_a_minus, delta_a_plus), units="dimensionless")
    # convert from log10(SFRD) to SFRD and carry the uncertainties
    SFR_minus = 10.0 ** (SFR + delta_SFR_minus)
    SFR_plus = 10.0 ** (SFR + delta_SFR_plus)
    SFR = 10.0 ** SFR
    SFR_scatter = unyt.unyt_array(
        (SFR - SFR_minus, SFR_plus - SFR), units="Msun/yr/Mpc**3"
    )
    SFR = unyt.unyt_array(SFR, units="Msun/yr/Mpc**3")

    processed.associate_x(
        a_bin,
        scatter=a_scatter,
        comoving=False,
        description="Cosmic scale factor",
    )
    processed.associate_y(
        SFR,
        scatter=SFR_scatter,
        comoving=False,
        description="Cosmic average star formation rate density",
    )

    z_minus = (z + delta_z_minus).min()
    z_plus = (z + delta_z_plus).max()
    processed.associate_redshift(0.5 * (z_minus + z_plus), z_minus, z_plus)
    processed.associate_plot_as(plot_as)

    output_path = f"{output_directory}/{output_filename}"

    if os.path.exists(output_path):
        os.remove(output_path)

    processed.write(filename=output_path)


# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Generate, format and save the Novak2017 data
cosmic_star_formation_history_novak()
