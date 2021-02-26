from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Chang2015.txt"
delimiter = None

output_filename = "Chang2015.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

processed = ObservationalData()
raw = np.loadtxt(input_filename, delimiter=delimiter)

comment = (
    "Assuming Chabrier IMF. z=0.01 - 0.2. No h-correction "
    "was required as data was supplied h-free. "
)
citation = "Chang et al. (2015) (SDSS)"
bibcode = "2013MNRAS.434..209B"
name = "Galaxy Stellar Mass - Galaxy Size from SDSS"
plot_as = "line"
redshift = 0.1
h = cosmology.h

log_M = raw.T[0]
M = unyt.unyt_array(10 ** (log_M), units=unyt.Solar_Mass)
SFR = unyt.unyt_array(10 ** raw.T[1], units=unyt.Solar_Mass / unyt.year)
sSFR = SFR / M

processed.associate_x(
    M, scatter=None, comoving=False, description="Galaxy Stellar Mass"
)
processed.associate_y(
    sSFR,
    scatter=None,
    comoving=False,
    description="Specific Star Formation Rate (sSFR)",
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
