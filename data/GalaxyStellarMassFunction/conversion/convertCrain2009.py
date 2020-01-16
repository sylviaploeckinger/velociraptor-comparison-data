from velociraptor.observations.objects import ObservationalData
from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Crain2009.txt"
delimiter = "\t"

output_filename = "Crain2009.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

processed = ObservationalData()
raw = np.loadtxt(input_filename, delimiter=delimiter)

comment = (
    f"Assuming Chabrier IMF. h-corrected for SWIFT using cosmology: {cosmology.name}."
)
citation = "Crain et al. (2009) (GIMIC)"
bibcode = "2009MNRAS.399.1773C"
name = "GSMF from GIMIC"
plot_as = "line"
redshift = 0.0
h = cosmology.h

log_M = raw.T[0]
M = 10 ** (log_M) * unyt.Solar_Mass / h
Phi = (10 ** raw.T[1] * (h ** 3)) * unyt.Mpc ** (-3)

processed.associate_x(M, scatter=None, comoving=True, description="Galaxy Stellar Mass")
processed.associate_y(Phi, scatter=None, comoving=True, description="Phi (GSMF)")
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
