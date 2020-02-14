from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Moustakas2013.txt"
delimiter = None

output_filename = "Moustakas2013.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

processed = ObservationalData()
raw = np.loadtxt(input_filename, delimiter=delimiter)

comment = (
    "Assuming Chabrier IMF. z=0.01 - 0.2. No h-correction "
    "was required as data was supplied h-free. "
    "Adopts cosmological parameters of h=0.7, omega0=0.3 "
    "omegaL=0.7."
)
citation = "Moustakas et al. (2013) (SDSS)"
bibcode = "2013ApJ...767...50M"
name = "Galaxy Stellar Mass - Passive Fraction from SDSS"
plot_as = "points"
redshift = 0.1
h = cosmology.h

log_M = raw.T[2]
M = unyt.unyt_array(10 ** (log_M), units=unyt.Solar_Mass)
passive_frac = unyt.unyt_array(raw.T[4] / raw.T[3], units="dimensionless")


processed.associate_x(
    M, scatter=None, comoving=False, description="Galaxy Stellar Mass"
)
processed.associate_y(
    passive_frac, scatter=None, comoving=False, description="Passive Fraction"
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
