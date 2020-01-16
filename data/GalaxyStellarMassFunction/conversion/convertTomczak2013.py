from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
# These lines are _required_ and you are required to use
# the cosmology specified (this is an astropy.cosmology
# instance)
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/AuthorYear.txt"
delimiter = None

output_filename = "AuthorYear.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

processed = ObservationalData()
raw = np.loadtxt(input_filename, delimiter=delimiter)

comment = f"Comment About Your Data. h-corrected for SWIFT using Cosmology: {cosmology.name}."
citation = "Author et al. (Year)"
bibcode = "Bibcode from ADS"
name = "Name of Plot"
plot_as = "points/line"
redshift = 0.000000
h = cosmology.h

x = convert_x_to_physical_units
y = convert_y_to_physical_units
# y_scatter should be a 1xN or 2xN array describing offsets from
# the median point 'y'
y_scatter = convert_y_scatter_to_physical_units

processed.associate_x(x, scatter=None, comoving=True, description="x Description")
processed.associate_y(y, scatter=y_scatter, comoving=True, description="y Description")
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
