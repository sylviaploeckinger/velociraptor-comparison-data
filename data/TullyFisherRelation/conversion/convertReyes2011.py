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

output_filename = "Reyes2011.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)


# Fits taken from eq. 39
stellar_masses = np.linspace(9.0, 11.0, 128)
v_max = 2.142 + 0.278 * (stellar_masses - 10.102)

stellar_masses = (
    unyt.unyt_array(10 ** stellar_masses, units=unyt.Solar_Mass)
    * kroupa_to_chabrier_mass
)
v_max = unyt.unyt_array(10 ** v_max, units=unyt.km / unyt.s)

# Meta-data
comment = (
    "Best-fit (eq. 39) to the data of 189 local galaxies. "
    "No cosmology correction needed as variables provided as physical. "
    f"Converted Kroupa to Chabrier IMF using ratio {kroupa_to_chabrier_mass}."
)
citation = "Reyes et al. (2011) (Fit)"
bibcode = "2011MNRAS.417.2347R"
name = "Fit to the stellar mass - v_80 (tully-fisher) relation at z=0."
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
