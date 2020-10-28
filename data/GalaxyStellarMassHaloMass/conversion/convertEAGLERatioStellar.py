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

# Redshifts at which to plot the data
redshifts = [0.0, 0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0]

haloes = ["M200", "MBN98"]
latex_names = ["$M_{200, {\\rm crit}}$", "$M_{\\rm BN98}$"]

# Create and save data for each halo definition and z in redshifts
for halo, latex_name in zip(haloes, latex_names):
    for z in redshifts:

        output_filename = f"Schaye2015_Ref_100_{f'{halo:s}'}_{f'z{z:07.3f}'.replace('.', 'p')}_RatioStellar.hdf5"
        output_directory = "../"

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
        name = f"Galaxy Stellar Mass / Halo mass - Stellar Mass relation at z={f'{z:3.1f}'}."
        plot_as = "line"
        redshift = z
        h = h_sim

        # Write everything
        processed = ObservationalData()
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
