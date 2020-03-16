from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Davies2016.txt"
delimiter = ","

output_filename = "Davies2016_z0p1.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

processed = ObservationalData()
raw = np.loadtxt(input_filename, delimiter=delimiter)

comment = (
    "Assuming Chabrier IMF (2003), and a standard LCDM cosmology with H=70km/s/Mpc "
    "and O_L = 0.7. Supplied h-free so no corrections have been made."
)
citation = "Davies et al. (2016)"
bibcode = "2016MNRAS.461..458D"
name = "Galaxy Stellar Mass-Star Formation Rate "
plot_as = "line"
redshift = 0.0
h_obs = 0.7
h = cosmology.h

M = raw.T[0] * unyt.Solar_Mass
SFR = raw.T[1] * unyt.Solar_Mass / unyt.year

processed.associate_x(M, scatter=None, comoving=True, description="Galaxy Stellar Mass")
processed.associate_y(
    SFR, scatter=None, comoving=True, description="Star Formation Rate"
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
