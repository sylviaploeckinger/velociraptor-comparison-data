from velociraptor.observations.objects import ObservationalData

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


def parse_latex_value(latex_string):
    """
    Take a LaTeX markup in the form ${value}_{-ve error}^{+ve error}$ and extract the
    numeric data from it.

    latex_string: The string to parse
    """
    values = re.findall("(?:-|\+)?(\d+.\d+|cdots)", latex_string)
    ret = []
    for v in values:
        # Missing data gets replaced with NaNs
        if "cdots" in v:
            return [np.nan, 0, 0]
        else:
            ret.append(float(v))
    return ret


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
    header_line_nos = [i for i, line in enumerate(lines) if "Tabulated" in line]
    header_line_nos.append(len(lines))

    # split the full list of lines into one block of lines per redshift bin
    split_lines = []
    for l1, l2 in pairwise(header_line_nos):
        split_lines.append(lines[l1:l2])

    # The values in the datafile are written in LaTeX markup
    converter_dict = dict(zip([1, 2], it.repeat(parse_latex_value)))

    # figure out the redshift bins
    z_bins_arr = np.zeros_like(split_lines)
    gsmf_arr = []
    for isl, lines in enumerate(split_lines):
        z_bins_arr[isl] = float(re.search("(\d.\d) <or=", lines[0]).group(1))
        # find the lines containing the actual data
        data_line_nos = [
            i for i, line in enumerate(lines) if line != "" and line[0].isdigit()
        ]
        data_slice = slice(data_line_nos[0], data_line_nos[-1])
        # The values in the datafile are written in val+-error tuples using LaTeX markup
        # We need to extract them manually
        data = np.zeros((data_slice.stop - data_slice.start, 4))
        for idt, line in enumerate(lines[data_slice]):
            line_split = line.split()
            extract = [float(line_split[0])]  # stellar mass bin
            extract.extend(parse_latex_value(line_split[1]))  # dusty GSMF
            extract.extend(parse_latex_value(line_split[2]))  # nondusty GSMF
            # The datafile tabulates the GSMF for dusty and non-dusty galaxies
            # We combine these to get a total GSMF, adding errors in quadrature
            tot_gsmf = extract[1] + extract[4]
            nve_err = tot_gsmf * np.sqrt(
                (extract[2] / extract[1]) ** 2 + (extract[5] / extract[4]) ** 2
            )
            pve_err = tot_gsmf * np.sqrt(
                (extract[3] / extract[1]) ** 2 + (extract[6] / extract[4]) ** 2
            )
            data[idt] = np.array([extract[0], tot_gsmf, nve_err, pve_err])
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

    comment = (
        "Assuming Chabrier IMF and Vmax selection, quoted redshift is lower bound of range. "
        f"h-corrected for SWIFT using Cosmology: {cosmology.name}."
    )
    citation = "Deshmukh et al. (2018)"
    bibcode = "2018ApJ...864..166D"
    name = "GSMF from SMUVS"
    plot_as = "points"
    redshift = z
    h = cosmology.h

    Mstar_bins = gsmf_and_Mstar_at_z[:, 0]
    M = 10 ** Mstar_bins * (h / ORIGINAL_H) ** (-2) * unyt.Solar_Mass
    # GSMF and errors are stored in datafile with units 10^-4 Mpc^-3 dex^-1
    Phi = 1e-4 * gsmf_and_Mstar_at_z[:, 1] * (h / ORIGINAL_H) ** 3 * unyt.Mpc ** (-3)
    # y_scatter should be a 1xN or 2xN array describing offsets from
    # the median point 'y'
    Phi_err = (
        1e-4 * gsmf_and_Mstar_at_z[:, 2:].T * (h / ORIGINAL_H) ** 3 * unyt.Mpc ** (-3)
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

input_filename = "../raw/Deshmukh2018.txt"

output_filename = "Deshmukh2018_{}.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# z_bins is a 1-D ndarray containing the lower edges of the redshift bins
# gsmf_and_Mstar is a list of 2D ndarrays, one per redshift
# Each contains four columns as follows:
# log(Mstar) bins, log(GSMF), GSMF +- errors
z_bins, gsmf_and_Mstar = load_file_and_split_by_z(input_filename)

for z, gsmf_and_Mstar_at_z in zip(z_bins, gsmf_and_Mstar):
    processed = process_for_redshift(z, gsmf_and_Mstar_at_z)

    output_path = f"{output_directory}/{output_filename.format(stringify_z(z))}"

    if os.path.exists(output_path):
        os.remove(output_path)

    processed.write(filename=output_path)
