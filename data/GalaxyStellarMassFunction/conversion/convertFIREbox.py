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

redshift = [0, 2, 4, 6, 8, 10]

comment = "Reweighted to account for cosmic variance in the halo mass function."
citation = f"Feldmann et al. (2022, FIREbox)"
bibcode = "2022arXiv220515325F"
name = f"Galaxy Stellar Mass Function (FIREbox)"
plot_as = "line"
multi_z = MultiRedshiftObservationalData()
multi_z.associate_citation(citation, bibcode)
multi_z.associate_name(name)
multi_z.associate_comment(comment)
multi_z.associate_cosmology(cosmology)

for z in redshift:
    input_filename = f"../raw/FIREbox_z{z:.0f}.txt"

    processed = ObservationalData()
    raw = np.loadtxt(input_filename)

    Mstar = 10.0 ** raw[:, 0]
    SMF = 10.0 ** raw[:, 1]
    Mstar = unyt.unyt_array(Mstar, units=unyt.Msun)
    SMF = unyt.unyt_array(SMF, units=1.0 / unyt.Mpc ** 3)

    processed.associate_x(
        Mstar, scatter=None, comoving=False, description="Galaxy Stellar Mass"
    )
    processed.associate_y(
        SMF,
        scatter=None,
        comoving=False,
        description="Galaxy Stellar Mass Function",
    )
    processed.associate_redshift(z)
    processed.associate_plot_as(plot_as)

    multi_z.associate_dataset(processed)

output_path = f"../FIREbox.hdf5"

if os.path.exists(output_path):
    os.remove(output_path)

multi_z.write(filename=output_path)
