from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys
import copy

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

output_filename_comw = "RemyRuyer2014_BPL_COMW.hdf5"
output_filename_coz = "RemyRuyer2014_BPL_COZ.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# The data is given in the form of a Schechter fit (Fig. 15 of the paper)
a = 2.21
alpha_hi = 1.00

# using X_CO,MW conversion factor
b_comw = 0.68
alpha_lo_comw = 3.08
oabundance_t_comw = 7.96

# using X_CO,MW conversion factor
b_coz = 0.96
alpha_lo_coz = 3.10
oabundance_t_coz = 8.10

# x value of the sun, 12+log10(O_solar/H_solar)
oabundance_sun = 8.69

# x-axis data, 12+log10(O/H)
oabundance = unyt.unyt_array(np.arange(7, 9.31, 0.01), units="dimensionless")

# Create the y-data for X_CO,MW case (BPL)
funcs = [
    lambda x: pow(10.0, -(a + alpha_hi * (oabundance_sun - x))),
    lambda x: pow(10.0, -(b_comw + alpha_lo_comw * (oabundance_sun - x))),
]
conds = [oabundance > oabundance_t_comw]

d2g_comw = unyt.unyt_array(
    np.piecewise(oabundance, conds, funcs), units="dimensionless"
)

# Create the y-data for X_CO,Z case (BPL)
funcs = [
    lambda x: pow(10.0, -(a + alpha_hi * (oabundance_sun - x))),
    lambda x: pow(10.0, -(b_coz + alpha_lo_coz * (oabundance_sun - x))),
]

conds = [oabundance > oabundance_t_coz]

d2g_coz = unyt.unyt_array(np.piecewise(oabundance, conds, funcs), units="dimensionless")

# Meta-data
comment = "Best-fit broken power law (BPL)"
citation_comw = "Rémy-Ruyer et al. [BPL, MW] (2014)"
citation_coz = "Rémy-Ruyer et al. [BPL, Z] (2014)"
bibcode = "2014A&A...563A..31R"
name_comw = "Dust-to-gas ratio as a function of metallicity assuming X_CO,MW conversion"
name_coz = "Dust-to-gas ratio as a function of metallicity assuming X_CO,Z conversion"
plot_as = "line"
redshift = 0.0
redshift_lower = 0.0
redshift_upper = 3.0
h = 0.7

# Write everything
outobj_comw = ObservationalData()
outobj_comw.associate_x(
    oabundance, scatter=None, comoving=True, description="12 + log10(O/H)"
)
outobj_comw.associate_y(
    d2g_comw,
    scatter=None,
    comoving=True,
    description="Dust-to-gas ratio (using X_CO,MW)",
)
outobj_comw.associate_citation(citation_comw, bibcode)
outobj_comw.associate_name(name_comw)
outobj_comw.associate_comment(comment)
outobj_comw.associate_redshift(redshift, redshift_lower, redshift_upper)
outobj_comw.associate_plot_as(plot_as)
outobj_comw.associate_cosmology(cosmology)

outobj_coz = copy.deepcopy(outobj_comw)
outobj_coz.associate_y(
    d2g_coz, scatter=None, comoving=True, description="Dust-to-gas ratio (using X_CO,Z)"
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
