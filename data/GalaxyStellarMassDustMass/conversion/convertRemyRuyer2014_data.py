from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys
import copy

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename_comw = "../raw/RemyRuyer2014_COMW.txt"
input_filename_coz = "../raw/RemyRuyer2014_COZ.txt"
delimiter = " "

output_filename_comw = "RemyRuyer2014_Data_COMW.hdf5"
output_filename_coz = "RemyRuyer2014_Data_COZ.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data, invert gas-to-dust to yield dust-to-gas
raw = np.loadtxt(input_filename_comw, delimiter=delimiter)
oabundance = raw[:, 0] * unyt.dimensionless  # 12 + log(O/H)
d2g_comw_med = pow(10, -raw[:, 1]) * unyt.dimensionless
d2g_comw_lo = pow(10, -raw[:, 2]) * unyt.dimensionless
d2g_comw_hi = pow(10, -raw[:, 3]) * unyt.dimensionless

raw = np.loadtxt(input_filename_coz, delimiter=delimiter)
d2g_coz_med = pow(10, -raw[:, 1]) * unyt.dimensionless
d2g_coz_lo = pow(10, -raw[:, 2]) * unyt.dimensionless
d2g_coz_hi = pow(10, -raw[:, 3]) * unyt.dimensionless

# Define the scatter as offset from the mean value
y_scatter_comw = unyt.unyt_array(
    (d2g_comw_med - d2g_comw_lo, d2g_comw_hi - d2g_comw_med)
)
y_scatter_coz = unyt.unyt_array((d2g_coz_med - d2g_coz_lo, d2g_coz_hi - d2g_coz_med))

# Meta-data
comment = "Median data binned by oxygen abundance 12+log10(O/H)"
citation_comw = "Rémy-Ruyer et al. [data, MW] (2014)"
citation_coz = "Rémy-Ruyer et al. [data, Z] (2014)"
bibcode = "2014A&A...563A..31R"
name_comw = "Dust-to-gas ratio as a function of metallicity assuming X_CO,MW conversion"
name_coz = "Dust-to-gas ratio as a function of metallicity assuming X_CO,Z conversion"
plot_as = "points"
redshift = 0.0
redshift_lower = 0.0
redshift_upper = 3.0
h = 0.7

# Write everything
outobj_comw = ObservationalData()
outobj_comw.associate_x(
    oabundance, scatter=None, comoving=True, description="Gas Phase 12 + log10(O/H)"
)
outobj_comw.associate_y(
    d2g_comw_med,
    scatter=y_scatter_comw,
    comoving=True,
    description="Dust-to-gas ratio (using X_CO,Z)",
)
outobj_comw.associate_citation(citation_comw, bibcode)
outobj_comw.associate_name(name_comw)
outobj_comw.associate_comment(comment)
outobj_comw.associate_redshift(redshift, redshift_lower, redshift_upper)
outobj_comw.associate_plot_as(plot_as)
outobj_comw.associate_cosmology(cosmology)

outobj_coz = copy.deepcopy(outobj_comw)
outobj_coz.associate_y(
    d2g_coz_med,
    scatter=y_scatter_coz,
    comoving=True,
    description="Dust-to-gas ratio (using X_CO,Z)",
)
outobj_comw.associate_citation(citation_coz, bibcode)
outobj_coz.associate_name(name_coz)

output_path_comw = f"{output_directory}/{output_filename_comw}"
output_path_coz = f"{output_directory}/{output_filename_coz}"

if os.path.exists(output_path_comw):
    os.remove(output_path_comw)

if os.path.exists(output_path_coz):
    os.remove(output_path_coz)

outobj_comw.write(filename=output_path_comw)
outobj_coz.write(filename=output_path_coz)
