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

input_filename = "../raw/Calette2021.txt"

output_template = "Calette2021_{}.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Meta-data
comment = (
    "Stellar Masses obtained assuming a Chabrier (2003) IMF. "
    "HI measurements from xGASS, following Catinella+18"
)

citation = "Calette et al. (2021)"
bibcode = "2021arXiv210401983C"
name = "Stellar mass - HI Gas to Stellar Mass ratio"
plot_as = "line"
redshift = 0.
h = h_sim

filetag = ['AllGalaxies', "Centrals", "Satellites"]
    
# Read the data
raw = np.loadtxt(input_filename)
M_star = pow(10., raw[:, 0]) * unyt.Solar_Mass

for i in range(3):
    MHI_per_Mstar =  pow(10., raw[:, i+1]) * unyt.dimensionless

    # Write everything
    processed = ObservationalData()
    processed.associate_x(
        M_star, scatter=None, comoving=True, description="Galaxy Stellar Mass"
    )
    processed.associate_y(
        MHI_per_Mstar, scatter=None, comoving=True, description="Stellar mass - HI Gas to Stellar Mass ratio"
    )
    processed.associate_citation(citation, bibcode)
    processed.associate_name(name)
    processed.associate_comment(comment)
    processed.associate_redshift(redshift, 0, 2)
    processed.associate_plot_as(plot_as)
    processed.associate_cosmology(cosmology)

    output_path = f"{output_directory}/{output_template.format(filetag[i])}"

    if os.path.exists(output_path):
        os.remove(output_path)

    processed.write(filename=output_path)
