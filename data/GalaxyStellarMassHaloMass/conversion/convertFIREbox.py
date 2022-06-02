from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

variations = {
    "Rvir": "Includes all stars within the virial radius.",
    "Rg": "Includes all stars within three times the half-light radius.",
}

for var in variations:
    input_filename = f"../raw/FIREbox_{var}.txt"

    processed = ObservationalData()
    raw = np.loadtxt(input_filename)

    comment = (
        f"Halos are identified using the BN98 virial overdensity. {variations[var]}"
    )
    citation = f"FIREbox ({var})"
    bibcode = "2022arXiv220515325F"
    name = f"Galaxy Stellar Mass Halo Mass relation ({var})"
    plot_as = "line"
    redshift = 0.0

    MBN98 = 10.0 ** raw[:, 0]
    Mratio = 10.0 ** raw[:, 1]
    Mstar = unyt.unyt_array(MBN98, units=unyt.Msun)
    Mratio = unyt.unyt_array(Mratio, units=unyt.dimensionless)

    processed.associate_x(
        Mstar, scatter=None, comoving=False, description="Galaxy Halo Mass (BN98)"
    )
    processed.associate_y(
        Mratio,
        scatter=None,
        comoving=False,
        description=f"Galaxy Stellar Mass to Halo Mass ratio ({var})",
    )
    processed.associate_citation(citation, bibcode)
    processed.associate_name(name)
    processed.associate_comment(comment)
    processed.associate_redshift(redshift)
    processed.associate_plot_as(plot_as)
    processed.associate_cosmology(cosmology)

    output_path = f"../FIREbox_{var}.hdf5"

    if os.path.exists(output_path):
        os.remove(output_path)

    processed.write(filename=output_path)
