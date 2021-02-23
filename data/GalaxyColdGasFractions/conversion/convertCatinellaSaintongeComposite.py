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

input_filename = "../raw/Catinella2018.txt"
input_filename_hi = "../raw/Catinella2018_HI.txt"
input_filename_h2 = "../raw/Saintonge2017_H2.txt"

output_filename = "CatinellaSaintongeComposite_{}.hdf5"
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

with open(input_filename_h2, "r") as f:
    lines = f.readlines()

#  load Saintonge+17 data
table_h2_index = [i for i, line in enumerate(lines) if "# x = " in line]

readrows = np.hstack([np.diff(table_h2_index).astype(int), None])

tables_h2 = []
for i in range(len(table_h2_index)):
    table_h2 = np.loadtxt(
        input_filename_h2, skiprows=table_h2_index[i], max_rows=readrows[i]
    )
    # set -99.9 values to NaNs
    table_h2[table_h2 == -99.9] = np.nan
    tables_h2.append(table_h2)

tables_h2 = [tables_h2[i] for i in [0, 1, 3, 2]]

# #  load Catinella+18 data
# table_hi_index = [i for i, line in enumerate(lines) if "# x = " in line]

# readrows = np.hstack([np.diff(table_hi_index).astype(int), None])

# tables_hi = []
# for i in range(len(table_hi_index)):
#     tables_hi.append(
#         np.loadtxt(input_filename, skiprows=table_hi_index[i], max_rows=readrows[i])
#     )


citation = "Catinella et al. (2018), z = 0"
comment = "HI+H2 vs stellar fractions at z=0, h-corrected for SWIFT using Cosmology: {cosmology.name}."
bibcode = "2018MNRAS.476..875C"
name = "H2/HI+H2 gas fractions from XGAS and XCOLDGAS"
plot_as = "points"
redshift = 0
redshift_lower = 0.0
redshift_upper = 4.0
h = cosmology.h

units = [
    pow(h / ORIGINAL_H, -2) * unyt.Solar_Mass,
    unyt.Solar_Mass * pow(unyt.kpc, -2),
    pow(unyt.yr, -1),
]

labels = [
    "Galaxy Stellar Mass",
    "Galaxy Central Stellar Surface density",
    "Galaxy Specific Star Formation Rate",
]

filetag = ["abcissa_M_star", "abcissa_mu_star", "abcissa_sSFR"]


for i in range(len(tables)):
    processed = ObservationalData()

    if i < 3:
        x_vals = 10 ** tables[i][:, 0] * units[i]
    else:
        x_vals = tables[i][:, 0] * units[i]

    # no x err
    x_err = x_vals * 0.0

    fneut = 10 ** tables[i][:, 1] * unitless
    log10_fh2 = np.interp(
        tables[i][:, 0],
        tables_h2[i][:, 0],
        tables_h2[i][:, 2],
        tables_h2[i][0, 2],
        tables_h2[i][-1, 2],
    )

    fh2 = 10 ** log10_fh2 * unitless

    fgas = fh2 / fneut
    fgas_err = fgas * 0.0

    processed.associate_x(x_vals, scatter=x_err, comoving=0, description=labels[i])
    processed.associate_y(
        fgas,
        scatter=fgas_err,
        comoving=0,
        description="Average galaxy cold gas to stellar fraction",
    )
    processed.associate_citation(citation, bibcode)
    processed.associate_name(name)
    processed.associate_comment(comment)
    processed.associate_redshift(redshift, redshift_lower, redshift_upper)
    processed.associate_plot_as(plot_as)
    processed.associate_cosmology(cosmology)

    output_path = f"{output_directory}/{output_filename.format(filetag[i])}"

    if os.path.exists(output_path):
        os.remove(output_path)

    processed.write(filename=output_path)
