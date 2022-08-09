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


# Fit directly from paper (Fig. 1).
alpha =   1.28
beta  = -11.40

log_sigma_H2  = np.linspace(np.log10(3.e6), 8.5)
log_sigma_sfr = alpha * log_sigma_H2 + beta

sigma_H2      = unyt.unyt_array(
    10 ** log_sigma_H2, units=unyt.Solar_Mass / (unyt.kpc * unyt.kpc)
)

sigma_sfr     = unyt.unyt_array(
    10 ** log_sigma_sfr, units=unyt.Solar_Mass / unyt.yr / (unyt.kpc * unyt.kpc)
)

# Meta-data
comment = (
    "Fit obtained from the paper, Equation 1 and Table 2 (see also Fig.1 )."
    "Median stellar mass of their sample is log Mstar = 10.35 and "
    "median pixel size is 0.3 kpc."
)
citation = "Abdurro'uf et al. (2022) (Fit)"
bibcode = "2019ApJ...876..155S"
name = "Fit to the resolved molecular Kennicutt-Schmidt relation at z=0."
plot_as = "line"
redshift = 0.0
h = h_sim   # They assume h = 0.7

# Write everything
processed = ObservationalData()
processed.associate_x(
    sigma_H2,
    scatter=None,
    comoving=False,
    description="$\\Sigma_{\\rm H_2}$",
)
processed.associate_y(
    sigma_sfr, scatter=None, comoving=False, description="$\\Sigma_{\\rm SFR}$"
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
