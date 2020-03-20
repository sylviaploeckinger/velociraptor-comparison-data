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

output_filename = "Sahu2019.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)


# Fit directly from paper.
log_velocity_dispersion = np.linspace(1.6, 2.7, 256)
velocity_dispersion = unyt.unyt_array(
    10 ** log_velocity_dispersion, units=unyt.km / unyt.s
)

log_black_hole_mass = 5.82 * (log_velocity_dispersion - np.log10(200)) + 8.17
black_hole_mass = unyt.unyt_array(10 ** log_black_hole_mass, units=unyt.Solar_Mass)


# Meta-data
comment = (
    "Fit obtained from the paper, Equation 3. Includes all galaxies "
    "from the sample of 137. Velocity dispersion was measured in fixed "
    "apertures of 0.595 h^-1 kpc using HYPERLEDA. Dispersions and "
    "masses provided h-free so no conversion applied."
)
citation = "Sahu et al. (2019) (Fit)"
bibcode = "2019ApJ...876..155S"
name = "Fit to the velocity-dispersion - black hole mass relation at z=0."
plot_as = "line"
redshift = 0.0
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    velocity_dispersion,
    scatter=None,
    comoving=False,
    description="Stellar velocity dispersion (~1kpc aperture)",
)
processed.associate_y(
    black_hole_mass, scatter=None, comoving=False, description="Black hole mass"
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
