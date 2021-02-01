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

output_filename = "Chartab2021_Data.hdf5"
output_directory = "../"
input_redshifts = [1.5, 2.3]

# Meta-data
comment = (
    "Fits to stacks obtained assuming a Chabrier IMF. "
    "The metallicity is expressed as 12 + log10(O/H). "
)
citation = "Chartab et al. (2021)"
bibcode = "https://arxiv.org/pdf/2101.01706.pdf"
name = "MOSDEF Survey: gas-phase metallicity of galaxies at 1.4 <= z <= 2.6"
plot_as = "line"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

multi_z = MultiRedshiftObservationalData()
multi_z.associate_citation(citation, bibcode)
multi_z.associate_name(name)
multi_z.associate_comment(comment)
multi_z.associate_cosmology(cosmology)
multi_z.associate_maximum_number_of_returns(1)

for filename, redshifts in zip(input_filenames, input_redshifts):
    processed = ObservationalData()
    raw = np.loadtxt(filename, delimiter=delimiter)

    plot_as = "line"
    redshift = redshifts
    h = cosmology.h

    log10_M_star_min = 9.5
    log10_M_star_max = 11.5

    if redshift == 1.5:

        Z_asm = 6.29
        alpha = 0.21
        M_star = np.arange(log10_M_star_min, log10_M_star_max, 0.2)
        Z_gas = (Z_asm + alpha * M_star) * unyt.dimensionless  # 12 + log(O/H)
        M_star = 10 ** M_star * unyt.Solar_Mass

    if redshift == 2.3:

        Z_asm = 4.84
        alpha = 0.35
        M_star = np.arange(log10_M_star_min, log10_M_star_max, 0.2)
        Z_gas = (Z_asm + alpha * M_star) * unyt.dimensionless  # 12 + log(O/H)
        M_star = 10 ** M_star * unyt.Solar_Mass

    processed.associate_x(
        M_star, scatter=None, comoving=True, description="Galaxy Stellar Mass"
    )
    processed.associate_y(
        Z_gas, scatter=None, comoving=True, description="Gas phase metallicity"
    )
    processed.associate_citation(citation, bibcode)
    processed.associate_name(name)
    processed.associate_comment(comment)
    processed.associate_redshift(redshift)
    processed.associate_plot_as(plot_as)
    processed.associate_cosmology(cosmology)

    multi_z.associate_dataset(processed)

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

multi_z.write(filename=output_path)
