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

output_filename = "Abdurrouf2022.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)


# Fit directly from paper (Fig. 10).
alpha = 1.65
beta = -11.81

log_sigma_gas = np.linspace(6.5, 8.5)
log_sigma_H2_over_sigma_HI = alpha * log_sigma_gas + beta

sigma_gas = unyt.unyt_array(
    10 ** log_sigma_gas, units=unyt.Solar_Mass / (unyt.kpc * unyt.kpc)
)

sigma_H2_over_sigma_HI = unyt.unyt_array(
    10 ** log_sigma_H2_over_sigma_HI, units=unyt.dimensionless
)

# Meta-data
comment = (
    "Fit obtained from the paper, Equation 1 and Fig.10 (middle panel)."
    "Median stellar mass of their sample is log Mstar = 10.35 and "
    "median pixel size is 0.3 kpc. The total gas surface density is "
    "Sigma_gas = 1.36 * Sigma_HI + Sigma_H2."
)
citation = "Abdurro'uf et al. (2022) (Fit)"
bibcode = "2022arXiv220708382A"
name = "Fit to the resolved HI to H2 transition at z=0."
plot_as = "line"
redshift = 0.0
h = h_sim  # They assume h = 0.7

# Write everything
processed = ObservationalData()
processed.associate_x(
    sigma_gas,
    scatter=None,
    comoving=False,
    description="$\\Sigma_{\\rm gas}$",
)
processed.associate_y(
    sigma_H2_over_sigma_HI,
    scatter=None,
    comoving=False,
    description="$\\Sigma_{\\rm H2} / \\Sigma_{\\rm HI}$",
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
