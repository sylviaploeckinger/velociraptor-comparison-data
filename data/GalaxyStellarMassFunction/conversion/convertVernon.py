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
name = "Galaxy Stellar Mass Function"
plot_as = "points"
redshift = 0.0
h_obs = 0.7
h = cosmology.h

log_mass=np.array([9.05,9.15,9.25,9.35,9.45,9.55,9.65,9.75,9.85,9.95,10.05,10.15,10.25,10.35,10.45,10.55,10.65,10.75,10.85,10.95,11.05,11.15,11.25,11.35,11.45,11.55,11.65,11.75,11.85,11.95,12.05,12.15])
log_smf=np.array([-2.032,-2.0518,-2.0698,-2.0948,-2.1342,-2.1678,-2.1977,-2.22,-2.2117,-2.227,-2.2371,-2.2419,-2.2524,-2.2671,-2.297,-2.3426,-2.4027,-2.4787,-2.5701,-2.6819,-2.8241,-2.9936,-3.2,-3.4466,-3.72,-4.0406,-4.3334,-4.6409,-4.9743,-5.505,-6.0078,-6.406])
log_std_dev=np.array([0.095047,0.085829,0.079139,0.072306,0.060862,0.060753,0.055938,0.050418,0.069255,0.06582,0.060716,0.056682,0.049267,0.043079,0.037501,0.035458,0.038148,0.043874,0.050008,0.061399,0.07296,0.086954,0.10556,0.13079,0.16787,0.24745,0.24816,0.27989,0.43343,0.46639,0.53379,0.41062])


mass = unyt.unyt_array(10**log_mass, units=unyt.Solar_Mass)
smf = unyt.unyt_array(10**log_smf, units=1 / unyt.Mpc**3)

scatter_down = 10**(log_smf) - 10**(log_smf - log_std_dev)
scatter_up = 10**(log_smf + log_std_dev) - 10**(log_smf)

smf_scatter = unyt.unyt_array(
    [scatter_down, scatter_up],
    units=1 / unyt.Mpc**3
)


processed.associate_x(
    mass, scatter=None, comoving=False, description="Galaxy Stellar Mass"
)
processed.associate_y(
    smf,
    scatter=smf_scatter,
    comoving=True,
    description="Galaxy Stellar Mass Function",
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
