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

    raw_file_name: the file name of the raw data file to extract the GSMF from
    """
    with open(raw_file_name, "r") as f:
        lines = f.readlines()

    z_bins_arr = []
    gsmf_arr = []

    for line in lines:
        lsplit = line.split("\t")  # file is tab-delimited
        # extract redshift info from header row
        if line.startswith("#"):
            for column_name in lsplit:
                z_match = re.search(r"z=([\d\.]+)", column_name)
                if z_match:
                    z_bins_arr.append(float(z_match.group(1)))
        else:
            # first field is the mass range
            mass_str = lsplit[0]
            mmin, mmax = map(float, mass_str.strip("[]").split("-"))

            # column starts with mass and mass bin width
            phi_vals = [0.5 * (mmin + mmax), 0.5 * (mmax - mmin)]

            # remaining fields are phi values, with +-error in brackets
            for phi_str in lsplit[1:]:
                phi_vals_match = re.search(
                    r"(-[\d\.]+)\((\+[\d\.]+) (-[\d\.]+)\)", phi_str
                )
                if phi_vals_match:
                    phi_vals.extend(map(float, phi_vals_match.groups()))
                else:  # data is missing
                    phi_vals.extend([np.nan] * 3)
            gsmf_arr.append(phi_vals)
    return np.array(z_bins_arr), np.array(gsmf_arr)


def process_for_redshift(z, gsmf_and_Mstar_at_z):
    """
    Output an HDF5 file containing the GSMF at a given redshift.

    z: the redshift to produce the GSMF for.
    gsmf_and_mstar_at_z: the array containing stellar mass bins and the GSMF at the
    chosen redshift
    """

    processed = ObservationalData()

    plot_as = "points"
    h = cosmology.h

    Mstar_bins = gsmf_and_Mstar_at_z[:, 0]
    M = 10 ** Mstar_bins * (h / ORIGINAL_H) ** (-2) * unyt.Solar_Mass

    M_err = (
        (10 ** Mstar_bins * np.log(10) * gsmf_and_Mstar_at_z[:, 1])
        * (h / ORIGINAL_H) ** (-2)
        * unyt.Solar_Mass
    )
    Phi = 10 ** gsmf_and_Mstar_at_z[:, 2] * (h / ORIGINAL_H) ** 3 * unyt.Mpc ** (-3)
    # y_scatter should be a 1xN or 2xN array describing offsets from
    # the median point 'y'
    # Errors are log error dz = 1/ln(10) dy/y
    # We want dy = y ln(10) dz
    Phi_err = (
        (
            10 ** gsmf_and_Mstar_at_z[:, 2][:, None]
            * np.abs(gsmf_and_Mstar_at_z[:, [4, 3]])
            * np.log(10)
        ).T
        * (h / ORIGINAL_H) ** 3
        * unyt.Mpc ** (-3)
    )

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

input_filename = "../raw/Gonzalez2011.txt"

output_filename = "Gonzalez2011.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

comment = (
    "Assuming Chabrier IMF (converted from Salpeter). "
    f"h-corrected for SWIFT using Cosmology: {cosmology.name}."
)
citation = "Gonzalez et. al. (2011)"
bibcode = "2011ApJ...735L..34G"
name = "GSMF from GOODS-S"

multi_z = MultiRedshiftObservationalData()
multi_z.associate_comment(comment)
multi_z.associate_name(name)
multi_z.associate_citation(citation, bibcode)
multi_z.associate_cosmology(cosmology)
multi_z.associate_maximum_number_of_returns(1)

# z_bins is a 1-D ndarray containing the lower edges of the redshift bins
# gsmf_and_Mstar is a list of 2D ndarrays, one per redshift
# Each contains five columns as follows:
# log(Mstar) bins, Mstar errors, log(GSMF), GSMF +- errors
z_bins, gsmf_and_Mstar = load_file_and_split_by_z(input_filename)

for i, z in enumerate(z_bins):
    gsmf_and_Mstar_at_z = np.concatenate(
        [gsmf_and_Mstar[:, :2], gsmf_and_Mstar[:, i * 3 + 2 : (i + 1) * 3 + 2]], axis=1
    )
    multi_z.associate_dataset(process_for_redshift(z, gsmf_and_Mstar_at_z))

multi_z.write(f"{output_directory}/{output_filename}")
