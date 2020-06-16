from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import re
import sys
import itertools as it

ORIGINAL_H = 0.7

unitless = (unyt.Solar_Mass / unyt.Solar_Mass).simplify()

# Exec the master cosmology file passed as first argument
# These lines are _required_ and you are required to use
# the cosmology specified (this is an astropy.cosmology
# instance)
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Catinella2018.txt"

output_filename = "Catinella2018_{}.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

with open(input_filename, "r") as f:
    lines = f.readlines()

# find header lines indicating the start of each block of dat
table_index = [i for i, line in enumerate(lines) if "# x = " in line]

readrows = np.hstack([np.diff(table_index).astype(int), None])

tables = []
for i in range(len(table_index)):
    tables.append(
        np.loadtxt(input_filename, skiprows=table_index[i], max_rows=readrows[i])
    )
    print(tables)

citation = "Catinella et al. (2018), z = 0"
comment = "HI to stellar fractions a z=0, h-corrected for SWIFT using Cosmology: {cosmology.name}."
bibcode = "2018MNRAS.476..875C"
name = "HI gas fractions from ALFALFA"
plot_as = "points"
redshift = 0
h = cosmology.h

units = [
    pow(h / ORIGINAL_H, -2) * unyt.Solar_Mass,
    unyt.Solar_Mass * pow(unyt.kpc, -2),
    pow(unyt.yr, -1),
    unitless,
]

labels = [
    "Galaxy Stellar Mass",
    "Galaxy Central Stellar Surface density",
    "Galaxy Specific Star Formation Rate",
    "Galaxy NUV-r colour",
]

filetag = ["abcissa_M_star", "abcissa_mu_star", "abcissa_sSFR", "abcissa_NUV_minus_r"]


for i in range(len(tables)):
    processed = ObservationalData()

    if i < 3:
        x_vals = 10 ** tables[i][:, 0] * units[i]
    else:
        x_vals = tables[i][:, 0] * units[i]

    # no x err
    x_err = x_vals * 0.0

    fhi = 10 ** tables[i][:, 1] * unitless

    fhi_plus_err = (10 ** (tables[i][:, 1] + tables[i][:, 2]) * unitless) - fhi
    fhi_minus_err = fhi - (10 ** (tables[i][:, 1] - tables[i][:, 2]) * unitless)
    fhi_err = np.row_stack([fhi_minus_err, fhi_plus_err])

    processed.associate_x(x_vals, scatter=x_err, comoving=0, description=labels[i])
    processed.associate_y(
        fhi,
        scatter=fhi_err,
        comoving=0,
        description="Average Galaxy HI to stellar fraction",
    )
    processed.associate_citation(citation, bibcode)
    processed.associate_name(name)
    processed.associate_comment(comment)
    processed.associate_redshift(redshift)
    processed.associate_plot_as(plot_as)
    processed.associate_cosmology(cosmology)

    output_path = f"{output_directory}/{output_filename.format(filetag[i])}"

    if os.path.exists(output_path):
        os.remove(output_path)

    processed.write(filename=output_path)
