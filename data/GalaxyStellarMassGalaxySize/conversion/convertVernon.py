"""
Conversion script for the parameter searching values.
"""

from velociraptor.observations.objects import ObservationalData
from velociraptor.tools.lines import binned_median_line

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

output_filename = "Vernon.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

processed = ObservationalData()


comment = (
    "Combination of various data sources, used for the calibration process "
    "in EAGLE-XL. Is h-free."
)
citation = "Calibration Data"
bibcode = "None available"
name = "Galaxy Stellar Mass-Galaxy Size"
plot_as = "points"
redshift = 0.0
h_obs = 0.7
h = cosmology.h

log_mass=np.array([7.1,7.3,7.5,7.7,7.9,8.1,8.3,8.5,8.7,8.9,9.1,9.3,9.5,9.7,9.9,10.1,10.3,10.5,10.7,10.9,11.1,11.3,11.5,11.7,11.9])
log_size=np.array([-0.28649,-0.22748,-0.16857,-0.10875,-0.05028,0.0082142,0.068261,0.13192,0.19398,0.2536,0.30157,0.3385,0.35724,0.36952,0.37548,0.38922,0.41625,0.46152,0.53019,0.61852,0.72904,0.8612,1.0152,1.1905,1.3851])
log_std_dev=np.array([0.48289,0.093203,0.083711,0.078168,0.070024,0.060474,0.049453,0.041253,0.033953,0.03754,0.037747,0.041652,0.041851,0.042617,0.040984,0.041501,0.04271,0.042675,0.040117,0.034548,0.02754,0.025737,0.046996,0.057768,0.10966])

mass = unyt.unyt_array(10**log_mass, units=unyt.Solar_Mass)
size = unyt.unyt_array(10**log_size, units=unyt.kpc)

scatter_down = 10**(log_size) - 10**(log_size - log_std_dev)
scatter_up = 10**(log_size + log_std_dev) - 10**(log_size)

size_scatter = unyt.unyt_array(
    [scatter_down, scatter_up],
    units=unyt.kpc
)


processed.associate_x(
    mass, scatter=None, comoving=False, description="Galaxy Stellar Mass"
)
processed.associate_y(
    size,
    scatter=size_scatter,
    comoving=False,
    description="Galaxy Half-Mass Radius",
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
