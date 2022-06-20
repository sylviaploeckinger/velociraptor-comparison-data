from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/FIREbox.txt"

processed = ObservationalData()
raw = np.loadtxt(input_filename)

comment = "Includes all galaxies in the simulation."
citation = "Feldmann et al. (2022, FIREbox)"
bibcode = "2022arXiv220515325F"
name = "Cosmic Average Stellar Mass Density (FIREbox)"
plot_as = "line"

a = raw[:, 0]
z = 1.0 / a - 1.0
rhostar = 10.0 ** raw[:, 1]
a = unyt.unyt_array(a, units=unyt.dimensionless)
rhostar = unyt.unyt_array(rhostar, units=unyt.Msun / unyt.Mpc ** 3)

processed.associate_x(
    a, scatter=None, comoving=False, description="Cosmic scale factor"
)
processed.associate_y(
    rhostar,
    scatter=None,
    comoving=False,
    description="Cosmic Stellar Mass Density",
)
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
zmin = z.min()
zmax = z.max()
processed.associate_redshift(0.5 * (zmin + zmax), zmin, zmax)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"../FIREbox.hdf5"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)
