from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Cosmology
h_obs = 0.7
h_sim = cosmology.h

output_filename = "Aird2015.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Data point taken from the paper conclusion
BHMD = unyt.unyt_array([4.20 * 1e5], "Msun / Mpc**3") / (h_obs / h_sim)
BHMD_m = unyt.unyt_array([0.14 * 1e5], "Msun / Mpc**3") / (h_obs / h_sim)
BHMD_p = unyt.unyt_array([0.29 * 1e5], "Msun / Mpc**3") / (h_obs / h_sim)

# Construct the error bar
y_scatter = unyt.unyt_array((BHMD_m, BHMD_p))

# Redshift of the data point
z = 0.01
a = unyt.unyt_array([1 / (1 + z)], "dimensionless")

# Meta-data
comment = (
    "Model fitted to X-ray luminosity data combining Chandra, ROSAT and ASCA surveys."
    f"h-corrected for SWIFT using cosmology: {cosmology.name}. "
)
citation = "Aird et al. (2015)"
bibcode = "2015MNRAS.451.1892A"
name = ""
plot_as = "points"
redshift = z
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(a, scatter=None, comoving=True, description="Scale-factor")
processed.associate_y(
    BHMD, scatter=y_scatter, comoving=True, description="Black-hole Mass Density"
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
