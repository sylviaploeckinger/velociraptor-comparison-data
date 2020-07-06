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

output_filename = "Andrews2013_Data.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Use the fit (equation 5) with the parameters from table 4
log10_M_star_min = 7.4
log10_M_star_max = 10.5
gamma = 0.640
log10_M_TO = 8.901
Z_asm = 8.798

M_TO = (10 ** log10_M_TO) * unyt.Solar_Mass

M_star = np.logspace(log10_M_star_min, log10_M_star_max, 100) * unyt.Solar_Mass
Z_star = (
    Z_asm - np.log10(1.0 + (M_TO / M_star) ** gamma)
) * unyt.dimensionless  # 12 + log(O/H)

# Convert the masses to Chabrier IMF
M_star = M_star * kroupa_to_chabrier_mass

# Meta-data
comment = (
    "Fits to the stacks obtained assuming a Kroupa IMF. "
    "The metallicity is expressed as 12 + log10(O/H). "
    "In these units the solar metallicity is 8.69."
    f"A correction of {kroupa_to_chabrier_mass} on the stellar "
    "masses has been applied to convert from Kroupa to the Chabrier IMF."
)
citation = "Andrews & Martini (2013) (SDSS)"
bibcode = "2013ApJ...765..140A"
name = "Stellar mass - Gas phase metallicity relation"
plot_as = "line"
redshift = 0.08
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star, scatter=None, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    Z_star, scatter=None, comoving=True, description="Gas phase metallicity"
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
