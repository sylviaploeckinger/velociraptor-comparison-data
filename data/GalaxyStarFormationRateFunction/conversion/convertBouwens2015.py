from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys
import re
import itertools as it


def pairwise(iterable):
    """
    return successive pairs of elements from the iterable
    i.e. (i0, i1), (i1, i2), (i2,i3), ...

    it: the iterable to consume.
    """
    a, b = it.tee(iterable)
    next(b, None)
    return zip(a, b)


def parse_latex_value(latex_string):
    """
    Take a LaTeX markup in the form ${value}_{-ve error}^{+ve error}$ and extract the
    numeric data from it.

    latex_string: The string to parse
    """
    values = re.findall("(\d+.\d+)", latex_string)
    ret = []
    for v in values:
        # Missing data gets replaced with NaNs
        ret.append(float(v))
    return ret


def load_file_and_split_by_z(raw_file_name):
    """
    Read the data file and do all the mucking around needed to extract a list of the
    redshift bins for which the SFRF is tabulated, along with the corresponding SFRF
    values and their errors.
    The number and spacing of the stellar mass bins vary with z; they are given in the
    first column of the returned array.

    raw_file_name: the file name of the raw data file to extract the SFRF from
    """
    with open(input_filename, "r") as f:
        lines = f.readlines()

    # find header lines indicating the start of each block of data
    header_line_nos = [i for i, line in enumerate(lines) if "z ∼" in line]
    header_line_nos.append(len(lines))

    # split the full list of lines into one block of lines per redshift bin
    split_lines = []
    for l1, l2 in pairwise(header_line_nos):
        split_lines.append(lines[l1:l2])

    z_bins_arr = [
        float(re.match("#z ∼ (\d.\d)", lines[lno]).group(1))
        for lno in header_line_nos[:-1]
    ]
    uv_lf_arr = []
    for isl, lines in enumerate(split_lines):
        uv_lf_arr.append(np.genfromtxt(lines, usecols=(0, 1, 3)))

    return z_bins_arr, uv_lf_arr


def process_for_redshift(z, sfrf_at_z):
    """
    Output an HDF5 file containing the SFRF at a given redshift.

    z: the redshift to produce the SFRF for.
    sfrf_at_z: the array containing SFR and Phi_SFR bins at the chosen redshift
    """

    processed = ObservationalData()

    comment = (
        "Assuming Chabrier IMF and Vmax selection. Includes dust corrections as "
        "described in Katsianis et. al. 2017, section 2."
        f"h-corrected for SWIFT using Cosmology: {cosmology.name}."
    )
    citation = "Bouwens et al. (2015)"
    bibcode = "2015ApJ...803...34B"
    name = "SFRF from HST Legacy Fields UV luminosity function"
    plot_as = "points"
    redshift = z
    h = cosmology.h

    SFR_bins = sfrf_at_z[:, 0]
    SFR = SFR_bins * unyt.Solar_Mass / unyt.year
    # SFRF and errors are stored in datafile with units 10^-2 Mpc^-3 dex^-1
    Phi = 1e-2 * sfrf_at_z[:, 1] * unyt.Mpc ** (-3)
    # y_scatter should be a 1xN or 2xN array describing offsets from
    # the median point 'y'
    Phi_err = 1e-2 * sfrf_at_z[:, 2].T * unyt.Mpc ** (-3)

    processed.associate_x(
        SFR, scatter=None, comoving=True, description="Star Formation Rate"
    )
    processed.associate_y(Phi, scatter=Phi_err, comoving=True, description="Phi (SFRF)")
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

input_filename = "../raw/Bouwens2015.txt"
delimiter = "\t"

output_filename = "Bouwens2015_{}.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

z_bins, UV_LF = load_file_and_split_by_z(input_filename)

for z, UV_LF_at_z in zip(z_bins, UV_LF):
    processed = process_for_redshift(z, UV_LF_at_z)

    output_path = f"{output_directory}/{output_filename.format(stringify_z(z))}"

    if os.path.exists(output_path):
        os.remove(output_path)

    processed.write(filename=output_path)
