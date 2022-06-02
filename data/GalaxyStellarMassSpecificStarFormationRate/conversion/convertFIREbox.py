from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

redshift = [0, 2]
variations = {
    "all": ("AllGalaxies", "Includes all galaxies, including quenched"),
    "SF": ("StarForming", "Only includes galaxies with a sSFR > 1.e-11 yr^-1"),
}

for z in redshift:
    for var in variations:
        input_filename = f"../raw/FIREbox_z{z:.0f}_{var}.txt"

        processed = ObservationalData()
        raw = np.loadtxt(input_filename)

        comment = variations[var][1]
        citation = f"FIREbox ({var})"
        bibcode = "2022arXiv220515325F"
        name = f"Galaxy Stellar Mass - Galaxy sSFR from FIREbox ({var})"
        plot_as = "points"
        redshift = z

        Mstar = 10.0 ** raw[:, 0]
        SFR = 10.0 ** raw[:, 1]
        sSFR = SFR / Mstar
        Mstar = unyt.unyt_array(Mstar, units=unyt.Msun)
        sSFR = unyt.unyt_array(sSFR, units=1.0 / unyt.yr)

        processed.associate_x(
            Mstar, scatter=None, comoving=False, description="Galaxy Stellar Mass"
        )
        processed.associate_y(
            sSFR,
            scatter=None,
            comoving=False,
            description="Specific Star Formation Rate (sSFR)",
        )
        processed.associate_citation(citation, bibcode)
        processed.associate_name(name)
        processed.associate_comment(comment)
        processed.associate_redshift(redshift)
        processed.associate_plot_as(plot_as)
        processed.associate_cosmology(cosmology)

        output_path = f"../FIREbox_z{z}_{variations[var][0]}.hdf5"

        if os.path.exists(output_path):
            os.remove(output_path)

        processed.write(filename=output_path)
