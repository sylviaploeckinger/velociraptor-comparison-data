from velociraptor.observations.objects import (
    ObservationalData,
    MultiRedshiftObservationalData,
)

import unyt
import numpy as np
import os
import sys


def cosmic_star_formation_history_gruppioni():

    # Meta-data
    name = f"Star formation rate density from Gruppioni et al. (2020)"
    comment = (
        "Uses the Chabrier initial mass function. "
        "Cosmology: H0=70, OmegaM=0.3, OmegaL=0.7"
    )

    citation = "Gruppioni et al. (2020)"
    bibcode = "2020A&A...643A...8G"
    plot_as = "points"
    output_filename = "Gruppioni2020.hdf5"
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
    data = np.loadtxt(f"../raw/sfr_gruppioni2020.dat")

    # Fetch the fields we need
    z_minus, z_plus = data[:, 0], data[:, 1]
    SFR, SFR_min, SFR_max = data[:, 2], data[:, 3], data[:, 4]

    z_bin = unyt.unyt_array(0.5 * (z_minus + z_plus), units="dimensionless")
    z_scatter = unyt.unyt_array(
        (z_bin - z_minus, z_plus - z_bin), units="dimensionless"
    )
    SFR_scatter = unyt.unyt_array(
        (SFR - SFR_min, SFR_max - SFR), units="Msun/yr/Mpc**3"
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

    z_minus = z_minus.min()
    z_plus = z_plus.max()
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
cosmic_star_formation_history_gruppioni()
