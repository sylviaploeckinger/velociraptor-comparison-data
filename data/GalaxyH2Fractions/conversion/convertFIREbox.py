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

comment = "Showing all central galaxies, using a 10 kpc aperture. The FIREbox results have been multiplied with 1.36 to take into account the contribution from Helium, as is done in observations."
citation = f"Feldmann et al. (2022, FIREbox)"
bibcode = "2022arXiv220515325F"
name = f"Galaxy Stellar Mass - H2 mass fraction (FIREbox)"
plot_as = "line"
redshift = 0.0

Mstar = 10.0 ** raw[:, 0]
MH2 = 10.0 ** raw[:, 1]
H2frac = MH2 / Mstar
Mstar = unyt.unyt_array(Mstar, units=unyt.Msun)
H2frac = unyt.unyt_array(H2frac, units=unyt.dimensionless)

processed.associate_x(
    Mstar, scatter=None, comoving=False, description="Galaxy Stellar Mass"
)
processed.associate_y(
    H2frac,
    scatter=None,
    comoving=False,
    description="Galaxy H2 mass fraction (10 kpc aperture)",
)
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(redshift)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"../FIREbox.hdf5"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)
