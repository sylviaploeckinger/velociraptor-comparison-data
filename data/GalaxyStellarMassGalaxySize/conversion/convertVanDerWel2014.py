from velociraptor.observations.objects import ObservationalData
from velociraptor.tools.lines import binned_median_line

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = {
    "red": "../raw/VanDerWel2014_red.txt",
    "blue": "../raw/VanDerWel2014_blue.txt",
}

output_filename = {"red": "VanDerWel_red.hdf5", "blue": "VanDerWel_blue.hdf5"}
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)


for color in ["red", "blue"]:

    processed = ObservationalData()

    comment = (
        "Assuming Chabrier IMF (2003), assumes a standard LCDM cosmology with "
        "H=71 km/s/Mpc and O_lambda=0.73."
    )
    citation = f"van der Wel et al. (2014, {color}-sample)"
    bibcode = "2014ApJ.788.28V"
    name = "Galaxy Stellar Mass-Galaxy Size"
    plot_as = "points"
    redshift = 0.25
    h_obs = 0.71
    h = cosmology.h

    raw = np.loadtxt(input_filename[color])

    M_star = 10 ** (raw[:, 0]) * unyt.Solar_Mass
    R_half_16 = 10 ** (raw[:, 1]) * unyt.kpc
    R_half_50 = 10 ** (raw[:, 2]) * unyt.kpc
    R_half_84 = 10 ** (raw[:, 3]) * unyt.kpc

    y_scatter = unyt.unyt_array((R_half_50 - R_half_16, R_half_84 - R_half_50))

    processed.associate_x(
        M_star, scatter=None, comoving=False, description="Galaxy Stellar Mass"
    )
    processed.associate_y(
        R_half_50,
        scatter=y_scatter,
        comoving=False,
        description="Galaxy Half-Mass Radius",
    )
    processed.associate_citation(citation, bibcode)
    processed.associate_name(name)
    processed.associate_comment(comment)
    processed.associate_redshift(redshift)
    processed.associate_plot_as(plot_as)
    processed.associate_cosmology(cosmology)

    output_path = f"{output_directory}/{output_filename[color]}"

    if os.path.exists(output_path):
        os.remove(output_path)

    processed.write(filename=output_path)
