from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import re
import sys


ORIGINAL_H = 0.70


def handle_bad_values(parts):
    """
    Convert "missing data" values to NaNs and fix the "errors" on these values getting
    misinterpreted as additional data points.

    parts: list of strings obtained by splitting a line from the datafile on "   "
    """
    n_parts = len(parts)
    to_remove = []
    for i in range(n_parts):
        if parts[i] == "-99":
            parts[i] = "NaN 0 0"
            to_remove.extend([i + 1, i + 2])
        elif parts[i] == "":
            to_remove.append(i)
    for idx in reversed(to_remove):
        del parts[idx]
    return parts


def load_file_and_split_by_z(raw_file_name):
    """
    Read the data file and do all the mucking around needed to extract lists of the
    redshift and stellar mass bins for which the GSMF is tabulated, along with the
    corresponding GSMF values and their errors.

    raw_file_name: the file name of the raw data file to extract the GSMF from
    """
    with open(raw_file_name, "r") as f:
        lines = f.readlines()
    # count lines that start with comment (#) or are blank
    n_header_lines = sum(re.match("#.+|^\s*$", l) is not None for l in lines)

    z_ranges = lines[4]
    # Each range of redshifts is separated by two or more spaces
    z_ranges = re.split("#?\s{2,}", z_ranges)[1:]
    z_bins_arr = np.asarray([float(z_rge.split()[0]) for z_rge in z_ranges])

    n_redshift_bins = len(z_bins_arr)
    n_stellar_mass_bins = len(lines) - n_header_lines
    gsmf_arr = np.zeros((n_redshift_bins, n_stellar_mass_bins, 3))

    mass_bins_arr = np.zeros(n_stellar_mass_bins)

    for ism, l in enumerate(lines[n_header_lines:]):
        # The GSMF values for each redshift bin are separated by two or more spaces
        parts = re.split("#?\s{2,}", l)
        # The first number on each line is the stellar mass bin
        mass_bins_arr[ism] = float(parts[0])

        if any(p == "-99" or p == "" for p in parts):
            # this indicates "bad value" and the errors will be given as "0"
            # but we will incorrectly register them as further data values
            # due to the way the file is formatted
            # so we need to sort this out first
            parts = handle_bad_values(parts)

        for iz, part in enumerate(parts[1:]):
            # each redshift bin has a 3-tuple of (GSMF, +ve error, -ve error)
            phi, errp, errn = map(float, part.split())
            gsmf_arr[iz, ism] = phi, errp, errn

    return z_bins_arr, mass_bins_arr, gsmf_arr


def process_for_redshift(z, mstar_bins, gsmf_at_z):
    """
    Output an HDF5 file containing the GSMF at a given redshift.

    z: the redshift to produce the GSMF for. The given value corresponds to the lower
    edge of a range in redshift, which has width 0.25 for z < 1.5 and 0.5 for z >= 1.5
    mstar_bins: the list of stellar mass bins for which the GSMF is tabulated
    gsmf_at_z: the slice of the GSMF array at the chosen redshift
    """

    processed = ObservationalData()

    comment = (
        "Assuming Chabrier IMF, quoted redshift is lower bound of range. "
        f"h-corrected for SWIFT using Cosmology: {cosmology.name}."
    )
    citation = "Tomczak et al (2013)"
    bibcode = "2014ApJ...783...85T"
    name = "GSMF from ZFOURGE/CANDELS"
    plot_as = "points"
    redshift = z
    h = cosmology.h

    M = 10 ** mstar_bins * (h / ORIGINAL_H) ** (-2) * unyt.Solar_Mass
    Phi = 10 ** gsmf_at_z[:, 0] * (h / ORIGINAL_H) ** 3 * unyt.Mpc ** (-3)
    # y_scatter should be a 1xN or 2xN array describing offsets from
    # the median point 'y'
    # Errors are log error dz = 1/ln(10) dy/y
    # We want dy = y ln(10) dz
    Phi_err = (
        (10 ** gsmf_at_z[:, 0][:, None] * np.log(10) * gsmf_at_z[:, [2, 1]]).T
        * (h / ORIGINAL_H) ** 3
        * unyt.Mpc ** (-3)
    )

    processed.associate_x(
        M, scatter=None, comoving=True, description="Galaxy Stellar Mass"
    )
    processed.associate_y(Phi, scatter=Phi_err, comoving=True, description="Phi (GSMF)")
    processed.associate_citation(citation, bibcode)
    processed.associate_name(name)
    processed.associate_comment(comment)
    processed.associate_redshift(redshift)
    processed.associate_plot_as(plot_as)
    processed.associate_cosmology(cosmology)

    return processed


def stringify_z(z):
    """
    Eagle-style text formatting of redshift label.
    Example: z=1.5 will be printed as z001p500.

    z: The redshift to produce a label for
    """
    whole = int(z)
    frac = int(1000 * (z - whole))
    return f"z{whole:03d}p{frac:03d}"


# Exec the master cosmology file passed as first argument
# These lines are _required_ and you are required to use
# the cosmology specified (this is an astropy.cosmology
# instance)
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Tomczak2013.txt"

output_filename = "Tomczak2013_{}.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# z_bins, mstar_bins are 1-D ndarrays containing the lower edges of the redshift bins,
# and the log(stellar mass) bins respectively
# gsmf is a 3-D ndarray with axes 0 and 1 corresponding to z and stellar mass bins
# Axis 2 ranges from 0..2 and contains log(GSMF), and the +- errors respectively
z_bins, mstar_bins, gsmf = load_file_and_split_by_z(input_filename)

for iz, z in enumerate(z_bins):
    processed = process_for_redshift(z, mstar_bins, gsmf[iz])

    output_path = f"{output_directory}/{output_filename.format(stringify_z(z))}"

    if os.path.exists(output_path):
        os.remove(output_path)

    processed.write(filename=output_path)
