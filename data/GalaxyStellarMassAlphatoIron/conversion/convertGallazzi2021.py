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
h_obs = 0.704  # WMAP7
Z_solar_obs = 0.02

input_filename = "../raw/Gallazzi_2021_ascii.txt"
delimiter = "\t"

output_filename = "Gallazzi2021_Data.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data
raw = np.loadtxt(input_filename, delimiter=delimiter)
M_star = 10 ** raw[:, 0] * unyt.Solar_Mass * (h_sim / h_obs) ** -2

# Correction factor due to the difference in (X_O/X_Fe)_Sun
# from Grevesse & Sauval (1993) to Asplund+ (2009)

O_over_H_Grevesse93 = 8.83  # Grevesse & Sauval
Fe_over_H_Grevesse93 = 7.5  # Grevesse & Sauval

O_over_H_Andres89 = 8.93
Fe_over_H_Andres89 = 7.51

O_over_H_Asplund09 = 8.69
Fe_over_H_Asplund09 = 7.50

O_over_Fe_solar_Grevesse93 = O_over_H_Grevesse93 - Fe_over_H_Grevesse93
O_over_Fe_solar_Andres89 = O_over_H_Andres89 - Fe_over_H_Andres89
O_over_Fe_solar_Asplund09 = O_over_H_Asplund09 - Fe_over_H_Asplund09

correction_Sun_O_over_Fe = O_over_Fe_solar_Grevesse93 - O_over_Fe_solar_Asplund09

Z_median = (raw[:, 1] + correction_Sun_O_over_Fe) * unyt.dimensionless
Z_lo = (raw[:, 2] + correction_Sun_O_over_Fe) * unyt.dimensionless
Z_hi = (raw[:, 3] + correction_Sun_O_over_Fe) * unyt.dimensionless

# Define the scatter as offset from the mean value
y_scatter = unyt.unyt_array((Z_median - Z_lo, Z_hi - Z_median))

# Meta-data
comment = (
    "Data obtained assuming a Chabrier IMF and h=0.7. "
    f"h-corrected for SWIFT using cosmology: {cosmology.name}. "
    "The metallicity is expressed as [alpha/Fe], using [alpha/Fe]solar =. "
    "The error bars given the 16th and 84th percentile of the distribution. "
    f"This has been corrected to use Z_solar={solar_metallicity} (Asplund+ 2009)"
)
citation = "Gallazzi et al. (2021)"
bibcode = "2021MNRAS.502...4457"
name = "Stellar mass - [alpha/Fe] relation"
plot_as = "line"
redshift = 0.1
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star, scatter=None, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    Z_median, scatter=y_scatter, comoving=True, description="Stellar metallicity"
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
