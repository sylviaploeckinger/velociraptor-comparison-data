from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import re
import sys
import itertools as it

ORIGINAL_H = 0.7

unitless = unyt.dimensionless

# Exec the master cosmology file passed as first argument
# These lines are _required_ and you are required to use
# the cosmology specified (this is an astropy.cosmology
# instance)
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Saintonge2017.txt"

output_filename = "Saintonge2017_{}.hdf5"
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
    table = np.loadtxt(input_filename, skiprows=table_index[i], max_rows=readrows[i])
    # set -99.9 values to NaNs
    table[table == -99.9] = np.nan
    tables.append(table)

citation = "Saintonge et al. (2017), z = 0"
comment = (
    "H2 gas fraction at z=0, h-corrected for SWIFT using Cosmology: {cosmology.name}."
)
bibcode = "2017ApJS..233...22S"
name = "H2 gas fractions from XCOLDGASS"
plot_as = "points"
redshift = 0
h = cosmology.h

units = [
    pow(h / ORIGINAL_H, -2) * unyt.Solar_Mass,
    unyt.Solar_Mass * pow(unyt.kpc, -2),
    unitless,
    pow(unyt.yr, -1),
]

labels = [
    "Galaxy Stellar Mass",
    "Galaxy Central Stellar Surface density",
    "Galaxy NUV-r colour",
    "Galaxy Specific Star Formation Rate",
]

filetag = ["abcissa_M_star", "abcissa_mu_star", "abcissa_NUV_minus_r", "abcissa_sSFR"]


for i in range(len(tables)):
    processed = ObservationalData()

    if i != 2:
        x_vals = 10 ** tables[i][:, 0] * units[i]
    else:
        x_vals = tables[i][:, 0] * units[i]

    fh2 = 10 ** tables[i][:, 2] * unitless

    fh2_plus_err = (10 ** (tables[i][:, 2] + tables[i][:, 3]) * unitless) - fh2
    fh2_minus_err = fh2 - (10 ** (tables[i][:, 2] - tables[i][:, 3]) * unitless)
    fh2_err = np.row_stack([fh2_minus_err, fh2_plus_err])

    processed.associate_x(x_vals, scatter=None, comoving=False, description=labels[i])
    processed.associate_y(
        fh2, scatter=fh2_err, comoving=False, description="Average Galaxy H2 fraction"
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
