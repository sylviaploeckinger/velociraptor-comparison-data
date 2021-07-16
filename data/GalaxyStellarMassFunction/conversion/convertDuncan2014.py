from velociraptor.observations.objects import (
    ObservationalData,
    MultiRedshiftObservationalData,
)

import unyt
import numpy as np
import os
import re
import sys
import itertools as it

ORIGINAL_H = 0.7


def pairwise(iterable):
    """
    return successive pairs of elements from the iterable
    i.e. (i0, i1), (i1, i2), (i2,i3), ...

    it: the iterable to consume.
    """
    a, b = it.tee(iterable)
    next(b, None)
    return zip(a, b)


def load_file_and_split_by_z(raw_file_name):
    """
    Read the data file and do all the mucking around needed to extract a list of the
    redshift bins for which the GSMF is tabulated, along with the corresponding GSMF
    values and their errors.
    The number and spacing of the stellar mass bins vary with z; they are given in the
    first column of the returned array.

    raw_file_name: the file name of the raw data file to extract the GSMF from
    """
    with open(raw_file_name, "r") as f:
        lines = f.readlines()

    # find header lines indicating the start of each block of data
    header_line_nos = [i for i, line in enumerate(lines) if "# z=" in line]
    header_line_nos.append(len(lines))

    # split the full list of lines into one block of lines per redshift bin
    split_lines = []
    for l1, l2 in pairwise(header_line_nos):
        split_lines.append(lines[l1:l2])

    # figure out the redshift bins
    z_bins_arr = []
    gsmf_arr = []
    for isl, lines in enumerate(split_lines):
        z_bins_arr.append(float(re.search("z=(\d)", lines[0]).group(1)))
        # find the lines containing the actual data
        data = np.loadtxt(lines)
        gsmf_arr.append(data)
    return z_bins_arr, gsmf_arr


def process_for_redshift(z, gsmf_and_Mstar_at_z):
    """
    Output an HDF5 file containing the GSMF at a given redshift.

    z: the redshift to produce the GSMF for. The given value corresponds to the lower
    edge of a range in redshift of width 0.5 for z < 3.0, and width 1.0 for z > 3.0
    gsmf_and_mstar_at_z: the array containing stellar mass bins and the GSMF at the
    chosen redshift
    """

    processed = ObservationalData()

    plot_as = "points"
    redshift = z
    h = cosmology.h

    Mstar_bins = 10 ** gsmf_and_Mstar_at_z[:, 0]
    M = Mstar_bins * (h / ORIGINAL_H) ** (-2) * unyt.Solar_Mass
    # Mass errors are log error dz = 1/ln(10) dy/y
    # We want dy = y ln(10) dz
    M_err = (
        Mstar_bins
        * np.log(10)
        * gsmf_and_Mstar_at_z[:, 1]
        * (h / ORIGINAL_H) ** (-2)
        * unyt.Solar_Mass
    )

    Phi = gsmf_and_Mstar_at_z[:, 2] * (h / ORIGINAL_H) ** 3 * unyt.Mpc ** (-3)
    # y_scatter should be a 1xN or 2xN array describing offsets from
    # the median point 'y'
    Phi_err = gsmf_and_Mstar_at_z[:, 3:].T * (h / ORIGINAL_H) ** 3 * unyt.Mpc ** (-3)

    processed.associate_x(
        M, scatter=M_err, comoving=True, description="Galaxy Stellar Mass"
    )
    processed.associate_y(Phi, scatter=Phi_err, comoving=True, description="Phi (GSMF)")
    processed.associate_citation(citation, bibcode)
    processed.associate_name(name)
    processed.associate_comment(comment)
    processed.associate_redshift(redshift)
    processed.associate_plot_as(plot_as)
    processed.associate_cosmology(cosmology)

    return processed


# Exec the master cosmology file passed as first argument
# These lines are _required_ and you are required to use
# the cosmology specified (this is an astropy.cosmology
# instance)
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Duncan2014.txt"

output_filename = "Duncan2014.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

comment = (
    "Assuming Chabrier IMF and Vmax selection."
    f"h-corrected for SWIFT using Cosmology: {cosmology.name}."
)
citation = "Duncan et al. (2014)"
bibcode = "2014MNRAS.444.2960D"
name = "GSMF from CANDELS/GOODS-S"

multi_z = MultiRedshiftObservationalData()
multi_z.associate_comment(comment)
multi_z.associate_name(name)
multi_z.associate_citation(citation, bibcode)
multi_z.associate_cosmology(cosmology)
multi_z.associate_maximum_number_of_returns(1)

# z_bins is a 1-D ndarray containing the redshift bins
# gsmf_and_Mstar is a list of 2D ndarrays, one per redshift
# Each contains five columns as follows:
# log(Mstar) bins, log(Mstar) error, log(GSMF), GSMF -+ errors

z_bins, gsmf_and_Mstar = load_file_and_split_by_z(input_filename)

for z, gsmf_and_Mstar_at_z in zip(z_bins, gsmf_and_Mstar):
    multi_z.associate_dataset(process_for_redshift(z, gsmf_and_Mstar_at_z))

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

multi_z.write(output_path)
