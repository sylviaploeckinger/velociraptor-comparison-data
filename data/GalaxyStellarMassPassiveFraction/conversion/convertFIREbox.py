from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = f"../raw/FIREbox_z0.txt"

processed = ObservationalData()
raw = np.loadtxt(input_filename)

comment = (
    "All galaxies are included, passive galaxies are defined as sSFR < 1.e-11 yr^-1"
)
citation = f"FIREbox"
bibcode = "2022arXiv220515325F"
name = f"Galaxy Stellar Mass - passive fraction (FIREbox)"
plot_as = "points"
redshift = 0.0

Mstar = 10.0 ** raw[:, 0]
pasfrac = raw[:, 1]
Mstar = unyt.unyt_array(Mstar, units=unyt.Msun)
pasfrac = unyt.unyt_array(pasfrac, units=unyt.dimensionless)

processed.associate_x(
    Mstar, scatter=None, comoving=False, description="Galaxy Stellar Mass"
)
processed.associate_y(
    pasfrac,
    scatter=None,
    comoving=False,
    description="Passive fraction",
)
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(redshift)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"../FIREbox_z0.hdf5"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)
