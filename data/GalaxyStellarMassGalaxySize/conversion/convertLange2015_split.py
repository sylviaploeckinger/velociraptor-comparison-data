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
    "red": "../raw/Lange2015HBand_red.txt",
    "blue": "../raw/Lange2015HBand_blue.txt",
}

output_filename = {"red": "Lange2015HBand_red.hdf5", "blue": "Lange2015HBand_blue.hdf5"}
output_directory = "../"

selection = {"red": "(u-r) > 1.5", "blue": "(u-r) < 1.5"}

if not os.path.exists(output_directory):
    os.mkdir(output_directory)


for color in ["red", "blue"]:

    processed = ObservationalData()

    comment = (
        "Assuming Chabrier IMF (2003), assumes a standard LCDM cosmology with "
        "H=70 km/s/Mpc and O_lambda=0.7."
    )
    citation = f"Lange et al. (2015, H-band, {selection[color]})"
    bibcode = "2015MNRAS.447.2603L "
    name = "Galaxy Stellar Mass-Galaxy Size"
    plot_as = "points"
    redshift = 0.05
    h_obs = 0.7
    h = cosmology.h

    raw = np.loadtxt(input_filename[color])

    M_star = raw[:, 0] * unyt.Solar_Mass
    R_half_16 = raw[:, 1] * unyt.kpc
    R_half_50 = raw[:, 2] * unyt.kpc
    R_half_84 = raw[:, 3] * unyt.kpc

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
