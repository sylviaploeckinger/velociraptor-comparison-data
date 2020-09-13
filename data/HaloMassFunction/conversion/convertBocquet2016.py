from velociraptor.observations.objects import ObservationalData
from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

output_filename = "Bocquet2016.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

data = np.loadtxt("../raw/HMF_Colossus.txt")

M_200 = unyt.unyt_array(data[:, 0], units=unyt.Solar_Mass)
Phi = unyt.unyt_array(data[:, 2], units=unyt.Mpc ** (-3))

comment = (
    "Halo mass functions at z=0 from the colossus package (Diemer+18) "
    "assuming a 'planck13' cosmology. "
    "Corrected to h-free units and dn/dlog10(M) (from dn/dln(M))."
)
citation = "Bocquet et al. (2016) (hydro)"
bibcode = "2016MNRAS.456.2361B"
name = "Halo mass function fit to the Magneticum simulations"
plot_as = "line"
redshift = 0.0
h = 0.6777

processed = ObservationalData()
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
