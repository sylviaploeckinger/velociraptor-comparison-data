from velociraptor.observations.objects import (
    ObservationalData,
    MultiRedshiftObservationalData,
)

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Cosmology
h_sim = cosmology.h

# Redshifts at which to plot the data
redshifts = [0.0, 0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0]
# Create the formatted version of the above array
redshift_header_info = ", ".join([f"{z:.1f}" for z in redshifts])

haloes = ["M200", "MBN98"]
latex_names = ["$M_{200, {\\rm crit}}$", "$M_{\\rm BN98}$"]

for halo, latex_name in zip(haloes, latex_names):

    # Meta-data
    comment = (
        "Medians obtained directly from the subfind catalogs. "
        "Central galaxies only. Halos with 0 stellar mass included. "
        f"Halo mass definition is {halo:s}. "
        "Stellar masses use a 30kpc spherical aperture. "
        "No cosmology correction needed."
    )
    citation = "EAGLE - L100N1504"
    bibcode = "2015MNRAS.446..521S"
    name = f"Galaxy Stellar Mass / Halo mass - Stellar Mass relation at z at z=[{redshift_header_info:s}]"
    plot_as = "line"
    h = h_sim

    # Store metadata at the top level

    multi_z = MultiRedshiftObservationalData()
    multi_z.associate_citation(citation, bibcode)
    multi_z.associate_name(name)
    multi_z.associate_comment(comment)
    multi_z.associate_cosmology(cosmology)
    multi_z.associate_maximum_number_of_returns(1)

    output_filename = f"Schaye2015_Ref_100_{f'{halo:s}'}_RatioStellar.hdf5"
    output_directory = "../"

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    for z in redshifts:

        # Create a single observational-data instance at redshift z
        processed = ObservationalData()

        if not os.path.exists(output_directory):
            os.mkdir(output_directory)

        data = np.loadtxt(
            f"../raw/EAGLE_SMHM_{f'{halo:s}'}_ratio_stellar_L100N1504_{f'z{z:03.1f}'.replace('.', 'p')}.txt"
        )

        M_star = unyt.unyt_array(data[:, 0], units=unyt.Solar_Mass)
        M_star_ratio_16 = unyt.unyt_array(data[:, 1], units=unyt.dimensionless)
        M_star_ratio_50 = unyt.unyt_array(data[:, 2], units=unyt.dimensionless)
        M_star_ratio_84 = unyt.unyt_array(data[:, 3], units=unyt.dimensionless)

        # Define the scatter as offset from the mean value
        y_scatter = unyt.unyt_array(
            (M_star_ratio_50 - M_star_ratio_16, M_star_ratio_84 - M_star_ratio_50)
        )

        # A fitting function gives us the data at any z. Hence, no need to have \Delta z
        redshift_lower, redshift_upper = [z, z]

        processed.associate_x(
            M_star,
            scatter=None,
            comoving=True,
            description=f"Galaxy Stellar Mass (30kpc, 3D)",
        )
        processed.associate_y(
            M_star_ratio_50,
            scatter=y_scatter,
            comoving=True,
            description="Galaxy Stellar Mass (30kpc, 3D) / Halo Mass ({latex_name:s})",
        )

        processed.associate_redshift(z, redshift_lower, redshift_upper)
        processed.associate_plot_as(plot_as)

        multi_z.associate_dataset(processed)

    output_path = f"{output_directory}/{output_filename}"

    if os.path.exists(output_path):
        os.remove(output_path)

    multi_z.write(filename=output_path)
