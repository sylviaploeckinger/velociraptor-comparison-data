from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys
import copy
import re
from velociraptor.tools.lines import binned_median_line

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Cosmologies
h_sim = cosmology.h
h_obs = h_sim

surveys = ["DGS", "KINGFISH"]

for surv in surveys:

    input_filename = f"../raw/Relano2020_{surv}.txt"

    output_directory = "../"
    output_filename = f"Relano2020_Data_{surv}.hdf5"

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    # Read the data (and convert ids to 0s for now)
    raw = np.genfromtxt(input_filename)
    M_star = pow(10, raw[:, 0]) * unyt.Solar_Mass
    S2L = pow(10, raw[:, 1]) / unyt.year

    lines = []
    labels = []
    inc = 0

    # Meta-data
    comment = (
        "Data obtained directly from Relano et al 2020, which uses",
        "A Chabrier et al (2003) IMF intrinsically.",
    )
    citation = f"{surv} (Relano et. al 2020)"
    bibcode = "10.1051/0004-6361/201937087"
    name = "Specific star formation rate - Grain Size Ratio Data"
    plot_as = "points"
    redshift = 0.02
    h = h_sim

    # Write everything
    processed = ObservationalData()
    processed.associate_x(
        M_star,
        scatter=None,
        comoving=True,
        description="Galaxy Specific Star Formation Rate",
    )
    processed.associate_y(
        S2L, scatter=None, comoving=True, description="Galaxy Dust Grain Size Ratio"
    )
    processed.associate_citation(citation, bibcode)
    processed.associate_name(name)
    processed.associate_comment(comment)
    processed.associate_redshift(redshift, 0, 0.5)
    processed.associate_plot_as(plot_as)
    processed.associate_cosmology(cosmology)

    output_path = f"{output_directory}/{output_filename}"

    if os.path.exists(output_path):
        os.remove(output_path)

    processed.write(filename=output_path)
