from velociraptor.observations.objects import (
    ObservationalData,
    MultiRedshiftObservationalData,
)

import unyt
import numpy as np
import os
import re
import sys


ORIGINAL_H = 0.6774
box_sizes = [50, 100, 300]


def load_file(file_name):
    """
    Read the data file and extract the redshift and stellar mass bins
    at which the GSMF is tabulated, along ith the corresponding values.

    file_name: the file name of the raw data to extract the GSMF from
    """
    with open(file_name, "r") as f:
        lines = f.readlines()

    header = list(filter(lambda l: l.startswith("#"), lines))
    nhead = len(header)

    # the last line of the header gives the redshift bins the GSMF is tabulated at
    z_str = re.findall("z=(\d.\d)", header[-1])
    z_bins = np.array([float(zs) for zs in z_str])

    raw_data = np.loadtxt(lines[nhead:])
    Mstar_bins = raw_data[:, 0]
    gsmf = raw_data[:, 1:].T  # transpose the GSMF so that first axis is redshift

    return z_bins, Mstar_bins, gsmf


def process_for_redshift(z, Mstar_bins, gsmf_at_z):
    """
    Output an `ObservationalData` instance containing the GSMF at a given redshift.

    z: the redshift to produce the GSMF for
    Mstar_bins: the list of stellar mass bins for which the GSMF is tabulated
    gsmf_at_z: the values of the GSMF at the chosen redshift
    """

    processed = ObservationalData()

    plot_as = "line"
    h = cosmology.h

    M = Mstar_bins * (h / ORIGINAL_H) ** (-2) * unyt.Solar_Mass
    M_err = None
    Phi = gsmf_at_z * (h / ORIGINAL_H) ** 3 * unyt.Mpc ** (-3)
    Phi_err = None

    processed.associate_x(
        M, scatter=M_err, comoving=True, description="Galaxy Stellar Mass"
    )
    processed.associate_y(Phi, scatter=Phi_err, comoving=True, description="Phi (GSMF)")
    processed.associate_redshift(z)
    processed.associate_plot_as(plot_as)

    return processed


# Exec the master cosmology file passed as first argument
# These lines are _required_ and you are required to use
# the cosmology specified (this is an astropy.cosmology
# instance)
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Pillepich2018_TNG{box}.txt"

output_filename = "Pillepich2018_TNG{box}.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

comment = (
    "Assuming Chabrier IMF, values calculated using 30 pkpc apertures. "
    "Obtained from TNG data portal. "
    f"h-corrected for SWIFT using Cosmology: {cosmology.name}."
)
citation = "Pillepich et al. (2018)"
bibcode = "2018MNRAS.475..648P"
name = "GSMF from Illustris-TNG{box}"

for box_sz in box_sizes:
    multi_z = MultiRedshiftObservationalData()
    multi_z.associate_comment(comment)
    multi_z.associate_name(name.format(box=box_sz))
    multi_z.associate_citation(citation, bibcode)
    multi_z.associate_cosmology(cosmology)
    multi_z.associate_maximum_number_of_returns(1)

    z_bins, Mstar_bins, gsmf = load_file(input_filename.format(box=box_sz))

    for z, gsmf_at_z in zip(z_bins, gsmf):
        multi_z.associate_dataset(process_for_redshift(z, Mstar_bins, gsmf_at_z))

    output_path = f"{output_directory}/{output_filename.format(box=box_sz)}"

    if os.path.exists(output_path):
        os.remove(output_path)

    multi_z.write(output_path)
