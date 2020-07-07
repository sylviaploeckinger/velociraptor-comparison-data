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

output_filename = "AvilaReese2008.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)


# Fit from table 1
log10_stellar_masses = np.linspace(8.8, 11.7, 128)
log10_v_max = -0.650 + 0.274 * log10_stellar_masses

# Apply IMF correction (Avila-Reese, private communication)
log10_stellar_masses = log10_stellar_masses - 0.15

# Apply additional correction to the masses from appendix of Li & White 2009
stellar_masses = 1e10 * 10.0 ** (((log10_stellar_masses - 10.0) - 0.130) / 0.922)

# Give proper units
stellar_masses = unyt.unyt_array(stellar_masses, units=unyt.Solar_Mass)
v_max = unyt.unyt_array(10.0 ** log10_v_max, units=unyt.km / unyt.s)

# Meta-data
comment = (
    "Fit obtained directly from paper using webplotdigitizer. "
    "No cosmology correction needed as variables provided as physical. "
    "Extracted using 76 galaxies. "
    "Converted Chabrier IMF and use an additional correction from "
    "Appendix A1 of Li & White (2009)."
)
citation = "Avila-Reese et al. (2008) (Fit)"
bibcode = "2008AJ....136.1340A"
name = "Fit to the stellar mass - vmax (tully-fisher) relation at z=0."
plot_as = "line"
redshift = 0.0
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    stellar_masses, scatter=None, comoving=False, description="Galaxy Stellar Mass"
)
processed.associate_y(v_max, scatter=None, comoving=False, description="Halo VMax")
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
