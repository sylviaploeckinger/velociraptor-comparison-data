from velociraptor.observations.objects import ObservationalData
from velociraptor.fitting_formulae.smhmr import moster_raw, behroozi_raw

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Cosmology
h_sim = cosmology.h

output_filename = "Moster2013RatioStellar.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)


halo_masses = np.logspace(4, 16, 512)
stellar_masses = moster_raw(0.0, halo_masses)

# Meta-data
comment = (
    "Fit obtained directly from paper using the smhmr module in "
    "velociraptor-python. No cosmology correction needed. "
    "Shows the ratio betweeen stellar mass and halo mass"
)
citation = "Moster et al. (2013)"
bibcode = "2013MNRAS.428.3121M"
name = "Fit to the stellar mass - stellar halo mass relation at z=0."
plot_as = "line"
redshift = 0.0
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    stellar_masses * unyt.Solar_Mass,
    scatter=None,
    comoving=False,
    description="Galaxy Stellar Mass",
)
processed.associate_y(
    (stellar_masses / halo_masses) * unyt.dimensionless,
    scatter=None,
    comoving=True,
    description="Galaxy Stellar Mass / Halo Mass",
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
