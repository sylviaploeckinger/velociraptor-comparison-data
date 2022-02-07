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

    # Load raw Enia2022 data
    data = np.loadtxt(f"../raw/sfr_novak2017.dat")

    # Fetch the fields we need
    z, z_minus, z_plus = data[:, 0], data[:, 1], data[:, 2]
    SFR, SFR_minus, SFR_plus = data[:, 3], data[:, 4], data[:, 5]

    z_bin = unyt.unyt_array(z, units="dimensionless")
    z_scatter = unyt.unyt_array((-z_minus, z_plus), units="dimensionless")
    # convert from log10(SFRD) to SFRD and carry the uncertainties
    SFR_minus = 10.0 ** (SFR + SFR_minus)
    SFR_plus = 10.0 ** (SFR + SFR_plus)
    SFR = 10.0 ** SFR
    SFR_scatter = unyt.unyt_array(
        (SFR - SFR_minus, SFR_plus - SFR), units="Msun/yr/Mpc**3"
    )
    SFR = unyt.unyt_array(SFR, units="Msun/yr/Mpc**3")

    processed.associate_x(
        z_bin,
        scatter=z_scatter,
        comoving=False,
        description="Cosmic redshift",
    )
    processed.associate_y(
        SFR,
        scatter=SFR_scatter,
        comoving=False,
        description="Cosmic average star formation rate density",
    )

    z_minus = (z + z_minus).min()
    z_plus = (z + z_plus).max()
    processed.associate_redshift(0.5 * (z_minus + z_plus), z_minus, z_plus)
    processed.associate_plot_as(plot_as)

    output_path = f"{output_directory}/{output_filename}"

    if os.path.exists(output_path):
        os.remove(output_path)

    processed.write(filename=output_path)


# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Generate, format and save the Enia2022 data
cosmic_star_formation_history_novak()
