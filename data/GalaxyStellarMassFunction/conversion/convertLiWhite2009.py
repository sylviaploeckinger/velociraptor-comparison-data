from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/LiWhite2009.txt"
delimiter = None

output_filename = "LiWhite2009.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

processed = ObservationalData()
raw = np.loadtxt(input_filename, delimiter=delimiter)

comment = (
    "Note stellar mass here is Msun/h^2. Converted. Poisson errors. "
    f"h-corrected for SWIFT using cosmology: {cosmology.name}."
)
citation = "Li and White (2009)"
bibcode = "2009MNRAS.398.2177L"
name = "GSMF from SDSS DR7"
plot_as = "points"
redshift = 0.07
h = cosmology.h

log_M = raw.T[0]
M = 10 ** (log_M) * unyt.Solar_Mass / (h ** 2)
# TODO: X errors
Phi = (10 ** raw.T[1]) * unyt.Mpc ** (-3) * h ** 3
Phi_scatter = None

processed.associate_x(M, scatter=None, comoving=True, description="Galaxy Stellar Mass")
processed.associate_y(Phi, scatter=Phi_scatter, comoving=True, description="Phi, GSMF")
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
