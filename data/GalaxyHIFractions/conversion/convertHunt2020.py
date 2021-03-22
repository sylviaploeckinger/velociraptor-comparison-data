from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Cosmologies
h_sim = cosmology.h

input_filename = "../raw/Hunt2020.txt"

output_filename = "Hunt2020_Data.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data
raw = np.loadtxt(input_filename)
M_star = pow(10.0, raw[:, 0]) * unyt.Solar_Mass
MHI_per_Mstar = pow(10.0, raw[:, 1]) * unyt.dimensionless
MHI_per_Mstar_lo = pow(10.0, raw[:, 2]) * unyt.dimensionless
MHI_per_Mstar_hi = pow(10.0, raw[:, 3]) * unyt.dimensionless

y_scatter = unyt.unyt_array(
    [MHI_per_Mstar - MHI_per_Mstar_lo, MHI_per_Mstar_hi - MHI_per_Mstar]
)

# Meta-data
comment = (
    "Stellar Masses obtained assuming a Chabrier (2003) IMF. "
    "local measurements decoupled from the Hubble flow (no h)."
    "HI measurements via 21cm emission in the MAGMA sample."
)

citation = "Hunt et al 2020 (MAGMA)"
bibcode = "2020A&A...643A.180H"
name = "Stellar mass - HI Gas to Stellar Mass ratio"
plot_as = "points"
redshift = 0.0
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star, scatter=None, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    MHI_per_Mstar,
    scatter=y_scatter,
    comoving=True,
    description="Stellar mass - HI Gas to Stellar Mass ratio",
)
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(redshift, 0, 2)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)
