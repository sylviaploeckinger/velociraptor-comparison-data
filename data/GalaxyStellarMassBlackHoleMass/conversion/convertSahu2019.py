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

input_filename = "../raw/Sahu2019.csv"
delimiter = ","

output_filenames = ["Sahu2019_ETG.hdf5", "Sahu2019_LTG.hdf5"]
output_names = ["ETG", "LTG"]
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data
raw_strings = np.loadtxt(
    input_filename, delimiter=delimiter, usecols=(1, 2), dtype=str, skiprows=1
).T
raw_floats = np.loadtxt(
    input_filename, delimiter=delimiter, usecols=(3, 4, 7, 8), dtype=float, skiprows=1
).T
exclusion, gal_type = raw_strings
log_mass_BH, err_log_mass_BH, log_mass_gal, err_log_mass_gal = raw_floats

valid_galaxies = exclusion == "i"

# Galaxy Stellar Masses
gal_stellar_mass = unyt.unyt_array(10 ** (log_mass_gal), units="Solar_Mass")
gal_stellar_mass_down = gal_stellar_mass - unyt.unyt_array(
    10 ** (log_mass_gal - err_log_mass_gal), units="Solar_Mass"
)
gal_stellar_mass_up = (
    unyt.unyt_array(10 ** (log_mass_gal + err_log_mass_gal), units="Solar_Mass")
    + gal_stellar_mass
)

# Black Hole Masses
BH_mass = unyt.unyt_array(10 ** (log_mass_BH), units="Solar_Mass")
BH_mass_down = BH_mass - unyt.unyt_array(
    10 ** (log_mass_BH - err_log_mass_BH), units="Solar_Mass"
)
BH_mass_up = (
    unyt.unyt_array(10 ** (log_mass_BH + err_log_mass_BH), units="Solar_Mass") + BH_mass
)

for output_name, output_filename in zip(output_names, output_filenames):
    mask = np.logical_and(valid_galaxies, gal_type == output_name)

    num_galaxies = np.sum(mask)
    # Meta-data
    comment = (
        "Data obtained assuming the total stellar mass is the same as the bulge mass. "
        f"No cosmology correction needed. This file includes only {output_name}s."
    )
    citation = f"Sahu (2019) ({output_name})"
    bibcode = "2019ApJ...887...10S"
    name = (
        "Black hole mass - stellar mass relation from "
        f"{num_galaxies} local {output_name}s."
    )
    plot_as = "points"
    # We purposely make this data show up not only a z=0 but also at higher z
    redshift_lower, redshift_upper = -0.1, 3.1
    redshift = 0.0
    h = h_sim

    # Need to mask arrays individually
    gal_stellar_mass_scatter = unyt.unyt_array(
        [gal_stellar_mass_down[mask], gal_stellar_mass_up[mask]]
    )
    BH_mass_scatter = unyt.unyt_array([BH_mass_down[mask], BH_mass_up[mask]])

    # Write everything
    processed = ObservationalData()
    processed.associate_x(
        gal_stellar_mass[mask],
        scatter=gal_stellar_mass_scatter,
        comoving=True,
        description="Galaxy Stellar Mass",
    )
    processed.associate_y(
        BH_mass[mask],
        scatter=BH_mass_scatter,
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
