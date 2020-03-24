from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Cosmologies
h_obs = 0.7
h_sim = cosmology.h

input_filename = "../raw/Tremonti_2004_ascii.txt"
delimiter = "\t"

output_filename = "Tremonti2004_Data.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data
raw = np.loadtxt(input_filename, delimiter=delimiter)
M_star = 10 ** raw[:, 0] * unyt.Solar_Mass * (h_sim / h_obs) ** -2
Z_median = raw[:, 3] * unyt.dimensionless  # 12 + log(O/H)
Z_lo = raw[:, 2] * unyt.dimensionless  # 12 + log(O/H)
Z_hi = raw[:, 4] * unyt.dimensionless  # 12 + log(O/H)

# Define the scatter as offset from the mean value
y_scatter = unyt.unyt_array((Z_median - Z_lo, Z_hi - Z_median))

# Meta-data
comment = (
    "Data obtained assuming a Kroupa IMF and h=0.7. "
    f"h-corrected for SWIFT using cosmology: {cosmology.name}. "
    "The metallicity is expressed as 12 + log10(O/H). "
    "In these units the solar metallicity is 8.69."
    "The error bars given the 16th and 84th percentile of the distribution."
)
citation = "Tremonti et al. (2004) (SDSS)"
bibcode = "2004ApJ...613..898T"
name = "Stellar mass - Gas phase metallicity relation"
plot_as = "line"
redshift = 0.1
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star, scatter=None, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    Z_median, scatter=y_scatter, comoving=True, description="Gas phase metallicity"
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
