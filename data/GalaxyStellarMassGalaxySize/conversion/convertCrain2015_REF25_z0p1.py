from velociraptor.observations.objects import ObservationalData
from velociraptor.tools.lines import binned_median_line

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Crain2015_REF25_z0p1.txt"
delimiter = " "

output_filename = "Crain2015_REF25_z0p1.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

processed = ObservationalData()
raw = np.loadtxt(input_filename, delimiter=delimiter)

comment = f"Assuming Chabrier IMF (2003), and EAGLE cosmology. h-corrected for SWIFT using cosmology: {cosmology.name}."
citation = "Crain et al. (2015) (EAGLE REF 25 Mpc)"
bibcode = "2015MNRAS.450.1937C"
name = "Galaxy Stellar Mass-Galaxy Size EAGLE NoAGN (25 Mpc)"
plot_as = "line"
redshift = 0.100639
h_obs = 0.7
h = cosmology.h

M = raw.T[0] * unyt.Solar_Mass * (h_obs / h)
R = raw.T[1] * unyt.kpc * (h_obs / h)

bins = unyt.unyt_array(np.logspace(8.5, 12, 20), units=unyt.Solar_Mass)

# Now bin the line
centers, median, deviation = binned_median_line(x=M, y=R, x_bins=bins)

processed.associate_x(
    centers, scatter=None, comoving=True, description="Galaxy Stellar Mass (30kpc, 3D)"
)
processed.associate_y(
    median,
    scatter=deviation,
    comoving=True,
    description="Galaxy Half-Mass Radius (30kpc 3D)",
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
