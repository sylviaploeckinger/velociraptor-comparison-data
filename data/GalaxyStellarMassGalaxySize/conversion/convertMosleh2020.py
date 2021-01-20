from velociraptor.observations.objects import (
    ObservationalData,
    MultiRedshiftObservationalData,
)
from velociraptor.tools.adaptive import create_adaptive_bins
from velociraptor.tools.lines import binned_median_line

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Mosleh2020.csv"
delimiter = ","

output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

raw = np.loadtxt(input_filename, delimiter=delimiter)

comment = (
    "Assuming Chabrier IMF (2003), exponentially declining SFH, and a standard"
    "LCDM cosmology with H=70km/s/Mpc "
    "and O_L = 0.7. Supplied h-free so no"
    "corrections have been made."
)
citation = "Mosleh et al. (2020)"
bibcode = "2020ApJ...905..170M"
name = "Galaxy Stellar Mass-Galaxy Size"
plot_as = "points"
h_obs = 0.7
h = cosmology.h

stellar_mass_bin_range = unyt.unyt_array([10 ** 9.6, 10 ** 11.6], "Solar_Mass")
number_of_bins = 10

z = raw.T[1]
M = raw.T[2] * unyt.Solar_Mass
R = raw.T[3] * unyt.kpc
e_R = raw.T[4] * unyt.kpc
sf = raw.T[5].astype(bool)

multi_z_sf = MultiRedshiftObservationalData()
multi_z_sf.associate_comment(f"{comment} Includes SFing galaxies only.")
multi_z_sf.associate_name(f"{name} (SF)")
multi_z_sf.associate_citation(f"{citation} (SF)", bibcode)
multi_z_sf.associate_cosmology(cosmology)

multi_z_nsf = MultiRedshiftObservationalData()
multi_z_nsf.associate_comment(f"{comment} Includes quiescent galaxies only.")
multi_z_nsf.associate_name(f"{name} (Q)")
multi_z_nsf.associate_citation(f"{citation} (Q)", bibcode)
multi_z_nsf.associate_cosmology(cosmology)

redshift_bins = [[0.3, 0.7], [0.7, 1.0], [1.0, 1.3], [1.3, 2.0]]

for redshift_bin in redshift_bins:
    processed_sf = ObservationalData()
    processed_nsf = ObservationalData()

    mask_sf = np.logical_and.reduce([z > redshift_bin[0], z < redshift_bin[1], sf])
    mask_nsf = np.logical_and.reduce([z > redshift_bin[0], z < redshift_bin[1], ~sf])

    stellar_mass_sf = M[mask_sf]
    stellar_mass_nsf = M[mask_nsf]

    _, bins_sf = create_adaptive_bins(
        stellar_mass_sf,
        lowest_value=stellar_mass_bin_range[0],
        highest_value=stellar_mass_bin_range[1],
        base_n_bins=number_of_bins,
        logarithmic=True,
        stretch_final_bin=True,
    )
    _, bins_nsf = create_adaptive_bins(
        stellar_mass_nsf,
        lowest_value=stellar_mass_bin_range[0],
        highest_value=stellar_mass_bin_range[1],
        base_n_bins=number_of_bins,
        logarithmic=True,
        stretch_final_bin=True,
    )

    centers_sf, medians_sf, deviations_sf = binned_median_line(
        x=stellar_mass_sf, y=R[mask_sf], x_bins=bins_sf, return_additional=False
    )

    centers_nsf, medians_nsf, deviations_nsf = binned_median_line(
        x=stellar_mass_nsf, y=R[mask_nsf], x_bins=bins_nsf, return_additional=False
    )

    processed_sf.associate_x(
        centers_sf, scatter=None, comoving=False, description="Galaxy Stellar Mass"
    )
    processed_sf.associate_y(
        medians_sf,
        scatter=deviations_sf,
        comoving=False,
        description="Galaxy Half-Mass Radius",
    )
    processed_sf.associate_redshift(0.5 * sum(redshift_bin), *redshift_bin)
    processed_sf.associate_plot_as(plot_as)

    processed_nsf.associate_x(
        centers_nsf, scatter=None, comoving=False, description="Galaxy Stellar Mass"
    )
    processed_nsf.associate_y(
        medians_nsf,
        scatter=deviations_nsf,
        comoving=False,
        description="Galaxy Half-Mass Radius",
    )
    processed_nsf.associate_redshift(0.5 * sum(redshift_bin), *redshift_bin)
    processed_nsf.associate_plot_as(plot_as)

    multi_z_sf.associate_dataset(processed_sf)
    multi_z_nsf.associate_dataset(processed_nsf)


output_path_sf = f"{output_directory}/Mosleh2020_SF.hdf5"
output_path_nsf = f"{output_directory}/Mosleh2020_Q.hdf5"

if os.path.exists(output_path_sf):
    os.remove(output_path_sf)

if os.path.exists(output_path_nsf):
    os.remove(output_path_nsf)

multi_z_sf.write(filename=output_path_sf)
multi_z_nsf.write(filename=output_path_nsf)
