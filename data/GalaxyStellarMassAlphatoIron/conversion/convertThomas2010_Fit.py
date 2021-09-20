from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Cosmology
h_sim = cosmology.h
output_filename = "Thomas2010.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Correction due to the difference in (X_Fe/X_O)_Sun
# From Anders and Grevesse (1989) to Asplund+ (2009)
O_over_Fe_solar_Andres89 = 0.717  # in log10
O_over_Fe_solar_Asplung09 = 0.64  # in log10
correction_Sun_O_over_Fe = O_over_Fe_solar_Andres89 - O_over_Fe_solar_Asplung09

# Conversion from log velocity dispersion to [alpha/Fe]
# (eq. 1 in Thomas et al. 2010)
func_alpha_over_Fe = (
    lambda log_sigma: -0.55 + 0.33 * log_sigma + correction_Sun_O_over_Fe
)

# Conversion from log Mstar to log velocity dispersion (eq. 2 in
# Thomas et al. 2005, https://arxiv.org/pdf/astro-ph/0410209)
func_log_sigma = lambda log_Mstar: (log_Mstar - 0.63) / 4.52

# X values of the data
M_star = np.logspace(np.log10(3e10), np.log10(3e11), 30)

# IMF correction (from M* Kroupa to M* Chabrier)
# (Lacey+ 2016, https://arxiv.org/pdf/1509.08473, table B2)
correction_IMF = 0.74 / 0.81

# Conversions
M_star_Kroupa_log = np.log10(M_star / correction_IMF)
sigma_log = func_log_sigma(log_Mstar=M_star_Kroupa_log)

# Y values of the data
alpha_over_Fe = func_alpha_over_Fe(log_sigma=sigma_log)

# Meta-data
comment = (
    "Fit to the the observed [alpha/Fe]∗ - M* relation for early-type galaxies. "
    "This fit is valid for a large range in galaxy masses (sigma > 100 km s-1 or "
    "Mdyn > 3e10 Msun) and all environmental densities. The original data is a "
    "magnitude-limited sample of 48 023 galaxies in the redshift range 0.05 <= "
    "z <= 0.1 with apparent r-band magnitude brighter than 16.8 from "
    "the SDSS Data Release 4. "
    "Solar abundances are converted from Anders and Grevesse (1989) to "
    "Asplund+ (2009)."
)
citation = "Thomas et al. 2010"
bibcode = "10.1111/j.1365-2966.2010.16427.x"
name = "Fit to [alpha/Fe]* vs. stellar mass relation at z=0.1"
plot_as = "line"
redshift_lower, redshift_upper = -0.1, 1.1
redshift = 0.1
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    unyt.unyt_array(M_star, "Msun"),
    scatter=None,
    comoving=False,
    description="Galaxy Stellar Mass",
)
processed.associate_y(
    unyt.unyt_array(alpha_over_Fe, "Dimensionless"),
    scatter=None,
    comoving=True,
    description="Black Hole Mass",
)
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(redshift, redshift_lower, redshift_upper)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)
