from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Bauer2013.txt"
delimiter = None

output_filename = "Bauer2013_StarForming.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

processed = ObservationalData()
raw = np.loadtxt(input_filename, delimiter=delimiter)

comment = (
    "Assuming Chabrier IMF. z=0.05 - 0.32. No h-correction "
    "was required as data was supplied h-free. "
    "This catalogue includes star-forming galaxies only."
)
citation = "Bauer et al. (2013) (GAMA, SF)"
bibcode = "2013MNRAS.434..209B"
name = "Galaxy Stellar Mass - Galaxy Size from GAMA (SF)"
plot_as = "points"
redshift = 0.2
h = cosmology.h

log_M = raw.T[0]
M = unyt.unyt_array(10 ** (log_M), units=unyt.Solar_Mass)
sSFR = unyt.unyt_array(10 ** raw.T[1], units=1 / unyt.year)
sSFR_scatter = unyt.unyt_array(
    [sSFR.value - 10 ** raw.T[3], 10 ** raw.T[2] - sSFR.value], units=1 / unyt.year
)

processed.associate_x(
    M, scatter=None, comoving=False, description="Galaxy Stellar Mass"
)
processed.associate_y(
    sSFR,
    scatter=sSFR_scatter,
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
