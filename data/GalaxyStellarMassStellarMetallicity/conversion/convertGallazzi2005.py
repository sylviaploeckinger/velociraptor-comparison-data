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
h_obs = 0.7

input_filename = "../raw/Gallazzi_2005_ascii.txt"
delimiter = "\t"

output_filename = "Gallazzi2005_Data.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data
raw = np.loadtxt(input_filename, delimiter=delimiter)
M_star = 10 ** raw[:, 0] * unyt.Solar_Mass * (h_sim / h_obs) ** -2
Z_median = 10 ** raw[:, 1] * unyt.dimensionless  # Z/Zsun
Z_lo = 10 ** raw[:, 2] * unyt.dimensionless  # Z/Zsun
Z_hi = 10 ** raw[:, 3] * unyt.dimensionless  # Z/Zsun

# Define the scatter as offset from the mean value
y_scatter = unyt.unyt_array((Z_median - Z_lo, Z_hi - Z_median))

# Meta-data
comment = (
    "Data obtained assuming a Kroupa IMF and h=0.7. "
    f"h-corrected for SWIFT using cosmology: {cosmology.name}. "
    "The metallicity is expressed in units of solar metallicity. "
    "The error bars given the 16th and 84th percentile of the distribution."
)
citation = "Gallazzi et al. (2005) (SDSS)"
bibcode = "2005MNRAS.362...41G"
name = "Stellar mass - Stellar metallicity relation"
plot_as = "line"
redshift = 0.1
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star, scatter=None, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    Z_median, scatter=y_scatter, comoving=True, description="Stellar metallicity"
)
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(redshift)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)
