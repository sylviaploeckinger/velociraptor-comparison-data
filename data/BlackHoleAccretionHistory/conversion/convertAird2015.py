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

input_filename = "../raw/Aird_2015.txt"

output_filename = "Aird2015.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data
raw = np.loadtxt(input_filename)
a = unyt.unyt_array(raw[:, 0], "dimensionless")
BHARD = unyt.unyt_array(raw[:, 1], "Msun / yr / Mpc**3")

# Convert to redshift
z = unyt.unyt_array(1 / a - 1.0, "dimensionless")

# Correct for cosmology
BHARD = BHARD * (h_sim / h_obs) ** -2

# Meta-data
comment = (
    "Model fitted to X-ray luminosity data combining Chandra, ROSAT and ASCA surveys."
    f"h-corrected for SWIFT using cosmology: {cosmology.name}. "
)
citation = "Aird et al. (2015)"
bibcode = "2015MNRAS.451.1892A"
name = "Redshift - Black-hole Mass Accretion Rate Density relation"
plot_as = "line"
redshift = z
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(a, scatter=None, comoving=True, description="Scale-factor")
processed.associate_y(
    BHARD, scatter=None, comoving=True, description="Black-hole Accretion Rate Density"
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
