from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Trujillo2020.txt"
delimiter = ","

output_filename = "Trujillo2020.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

processed = ObservationalData()
raw = np.loadtxt(input_filename, delimiter=delimiter)

comment = f"Assuming Chabrier IMF (2003), and a standard LCDM cosmology with H=70km/s/Mpc and O_L = 0.7. h-corrected for SWIFT using cosmology: {cosmology.name}."
citation = "Trujillo et al. (2020)"
bibcode = "2020arXiv200102689T"
name = "Galaxy Stellar Mass-Galaxy Size"
plot_as = "points"
redshift = 0.0
h_obs = 0.7
h = cosmology.h

M = raw.T[0] * unyt.Solar_Mass
R = raw.T[1] * unyt.kpc

processed.associate_x(M, scatter=None, comoving=True, description="Galaxy Stellar Mass")
processed.associate_y(
    R, scatter=None, comoving=True, description="Galaxy Half-Mass Radius"
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
