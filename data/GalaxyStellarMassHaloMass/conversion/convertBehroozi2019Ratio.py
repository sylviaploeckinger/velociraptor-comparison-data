from velociraptor.observations.objects import ObservationalData
from velociraptor.fitting_formulae.smhmr import behroozi_2019_raw

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Cosmology
h_sim = cosmology.h

output_filename = "Behroozi2019Ratio.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Redshift at which to plot the data
redshift = 0.0

M_vir = np.logspace(7, 15, 512)

# Note that the function below actually takes M_{vir, peak}, and we are
# ignoring this fact
M_star = behroozi_2019_raw(redshift, M_vir)

# Correct M_vir --> M_200 (fit from EAGLE cosmology - cosmology dependence is very weak)
M_200 = M_vir / 1.2

# Meta-data
comment = (
    "The data is taken from https://www.peterbehroozi.com/data.html"
    "Median Fit to the raw data including both satellites and centrals."
    "The stellar mass is the true stellar mass (i.e. w/o observational corrections)"
    "The halo mass is the peak halo mass that follows the Bryan & Norman (1998)"
    "spherical overdensity definition"
    "Cosmology: Omega_m=0.307, Omega_lambda=0.693, h=0.678, sigma_8=0.823, n_s=0.96"
    "Shows the ratio between stellar mass and halo mass as a function of halo mass"
)
citation = "Behroozi et al. (2019)"
bibcode = "2019MNRAS.488.3143B"
name = "Fit to the stellar mass - halo mass relation at z={:.1}".format(redshift)
plot_as = "line"
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_200 * unyt.Solar_Mass,
    scatter=None,
    comoving=True,
    description="Halo Mass ($M_{200, {\rm crit}}$)",
)
processed.associate_y(
    (M_star / M_200) * unyt.dimensionless,
    scatter=None,
    comoving=True,
    description="Galaxy Stellar Mass / Halo Mass ($M_* / M_{200, {\rm crit}}$)",
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
