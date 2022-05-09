from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Cosmology
h_sim = cosmology.h

output_filename_distant = "Popping2022.hdf5"
output_filename_local = "DeVis2019.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)


# Fit from table 1
oabundance = np.linspace(6, 10, 128)
log10_d2g_distant = 1.3 * oabundance - 13.72
log10_d2g_local = 2.45 * oabundance - 23.30
d2g_distant = pow(10., log10_d2g_distant)
d2g_local = pow(10., log10_d2g_local)

# Give proper units
oabundance = unyt.unyt_array(oabundance, units=None)
d2g_distant = unyt.unyt_array(d2g_distant, units=None)
d2g_local = unyt.unyt_array(d2g_local, units=None)

# Meta-data
comment = (
    "Fit obtained directly from paper Table 1"
)
citation_distant = "Popping & Peroux (2022) (Fit)"
citation_local = "De Vis et al. (2019) (Fit)"
bibcode = "arXiv:2203.03686"
name_distant = "Fit to the [O/H] - Dust-to-gas ratio relation at z>0.5"
name_local = "Fit to the [O/H] - Dust-to-gas ratio relation at z=0"
plot_as = "line"
redshift_distant = 2.
redshift_local = 0.
h = h_sim

# Write distant
processed_distant = ObservationalData()
processed_distant.associate_x(
    oabundance, scatter=None, comoving=False, description="[O/H]"
)
processed_distant.associate_y(d2g_distant, scatter=None, comoving=False, description="Dust-to-Gas Ratio")
processed_distant.associate_citation(citation_distant, bibcode)
processed_distant.associate_name(name_distant)
processed_distant.associate_comment(comment)
processed_distant.associate_redshift(redshift_distant, 0.5, 5.5)
processed_distant.associate_plot_as(plot_as)
processed_distant.associate_cosmology(cosmology)

output_path = f"{output_directory}/{output_filename_distant}"

if os.path.exists(output_path):
    os.remove(output_path)

processed_distant.write(filename=output_path)

# Write local
processed_local = ObservationalData()
processed_local.associate_x(
    oabundance, scatter=None, comoving=False, description="[O/H]"
)
processed_local.associate_y(d2g_local, scatter=None, comoving=False, description="Dust-to-Gas Ratio")
processed_local.associate_citation(citation_local, bibcode)
processed_local.associate_name(name_local)
processed_local.associate_comment(comment)
processed_local.associate_redshift(redshift_local, 0., 0.5)
processed_local.associate_plot_as(plot_as)
processed_local.associate_cosmology(cosmology)

output_path = f"{output_directory}/{output_filename_local}"

if os.path.exists(output_path):
    os.remove(output_path)

processed_local.write(filename=output_path)
