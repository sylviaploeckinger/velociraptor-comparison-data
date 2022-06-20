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

redshift = [0, 2]
variations = {
    "all": ("AllGalaxies", "Includes all galaxies, including quenched"),
    "SF": ("StarForming", "Only includes galaxies with a sSFR > 1.e-11 yr^-1"),
}

for var in variations:

    comment = variations[var][1]
    citation = f"Feldmann et al. (2022, FIREbox, {var})"
    bibcode = "2022arXiv220515325F"
    name = f"Galaxy Stellar Mass - Galaxy sSFR from FIREbox ({var})"
    plot_as = "line"

    multi_z = MultiRedshiftObservationalData()
    multi_z.associate_citation(citation, bibcode)
    multi_z.associate_name(name)
    multi_z.associate_comment(comment)
    multi_z.associate_cosmology(cosmology)

    for z in redshift:
        input_filename = f"../raw/FIREbox_z{z:.0f}_{var}.txt"

        processed = ObservationalData()
        raw = np.loadtxt(input_filename)

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
        processed.associate_redshift(z)
        processed.associate_plot_as(plot_as)

        multi_z.associate_dataset(processed)

    output_path = f"../FIREbox_{variations[var][0]}.hdf5"

    if os.path.exists(output_path):
        os.remove(output_path)

    multi_z.write(filename=output_path)
