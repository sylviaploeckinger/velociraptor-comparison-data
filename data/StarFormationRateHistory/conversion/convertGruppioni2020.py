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

    # Load raw Gruppioni2020 data
    data = np.loadtxt(f"../raw/sfr_gruppioni2020.dat")

    # Fetch the fields we need
    z_minus, z_plus = data[:, 0], data[:, 1]
    SFR, SFR_min, SFR_max = data[:, 2], data[:, 3], data[:, 4]

    z = 0.5 * (z_minus + z_plus)

    a = 1.0 / (1.0 + z)
    a_minus = 1.0 / (1.0 + z_plus)
    a_plus = 1.0 / (1.0 + z_minus)

    a_bin = unyt.unyt_array(a, units="dimensionless")
    a_scatter = unyt.unyt_array((a - a_minus, a_plus - a), units="dimensionless")
    SFR_scatter = unyt.unyt_array(
        (SFR - SFR_min, SFR_max - SFR), units="Msun/yr/Mpc**3"
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

# Generate, format and save the Gruppioni2020 data
cosmic_star_formation_history_gruppioni()
