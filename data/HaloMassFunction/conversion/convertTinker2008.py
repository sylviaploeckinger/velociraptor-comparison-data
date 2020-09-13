from velociraptor.observations.objects import ObservationalData
from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/HMF_Colossus.txt"
delimiter = "\t"

output_filename = "Tinker2008.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

processed = ObservationalData()
raw = np.loadtxt(input_filename, delimiter=delimiter)

M_200 = raw[:, 0]
Phi = raw[:, 1]

comment = (
    "Halo mass functions at z=0 from the colossus package (Diemer+18) assuming a 'planck13' cosmology. "
    "Corrected to h-free units and dn/dlog10(M) (from dn/dln(M))."
)
citation = "Tinker et al. (2008)"
bibcode = "2008ApJ...688..709T"
name = "Halo mass function fit to simulations"
plot_as = "line"
redshift = 0.0
h = 0.6777

processed.associate_x(M_200, scatter=None, comoving=True, description="Halo Mass")
processed.associate_y(Phi, scatter=None, comoving=True, description="Phi (HMF)")
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
