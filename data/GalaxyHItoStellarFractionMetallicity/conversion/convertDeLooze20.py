from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import re
import sys
from velociraptor.tools.lines import binned_median_line

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/JINGLE_MHIMstar_vs_Metallicity.dat"
delimiter = " "

output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

lines = []
labels = []
inc = 0

# find all the sub headers in the file (tthey contain a ':')
for inc, line in enumerate(open(input_filename)):
    if re.findall(r"\:", line):
        # store sub-head lines and survey names
        lines.append(inc)
        labels.append(line.split(":")[0])
    inc += 1
lines.append(inc)

x_all = []
y_all = []

# Loop over different survey data sets
for i in range(len(labels)):
    # read each surveys data and write separate files
    raw = np.genfromtxt(
        input_filename,
        skip_header=lines[i] + 1,
        skip_footer=inc - lines[i + 1],
        dtype=float,
        usecols=[1, 2, 3, 4, 5, 6],
    )
    label = "_".join(labels[i].split(" "))
    output_filename = f"DeLooze20_individual_{label}.hdf5"

    y_all.append(pow(10, raw[:, 0]))
    x_all.append(raw[:, 1])

    logMHIMstar_med = pow(10, raw[:, 0]) * unyt.dimensionless
    logMHIMstar_lo = pow(10, raw[:, 2]) * unyt.dimensionless
    logMHIMstar_hi = pow(10, raw[:, 4]) * unyt.dimensionless
    oabundance_med = raw[:, 1] * unyt.dimensionless
    oabundance_lo = raw[:, 3] * unyt.dimensionless
    oabundance_hi = raw[:, 5] * unyt.dimensionless

    # Define the scatter as offset from the mean value
    y_scatter = unyt.unyt_array(
        (logMHIMstar_med - logMHIMstar_lo, logMHIMstar_hi - logMHIMstar_med)
    )
    x_scatter = unyt.unyt_array(
        (oabundance_med - oabundance_lo, oabundance_hi - oabundance_med)
    )

    # Meta-data
    comment = f"values for individual galaxies derived from {label} survey"
    citation = f"{label} compiled by De Looze et al. (2020)"
    bibcode = "2020MNRAS.496.3668D"
    name = "MHI/Mstar as a function of 12 + log10(O/H)"
    plot_as = "points"
    redshift = 0.0
    redshift_lower = 0.0
    redshift_upper = 3.0
    h = 0.7

    # Write everything
    outobj = ObservationalData()
    outobj.associate_x(
        oabundance_med,
        scatter=x_scatter,
        comoving=True,
        description="Gas phase 12 + log10(O/H)",
    )
    outobj.associate_y(
        logMHIMstar_med,
        scatter=y_scatter,
        comoving=True,
        description="HI mass to stellar mass ratio",
    )
    outobj.associate_citation(citation, bibcode)
    outobj.associate_name(name)
    outobj.associate_comment(comment)
    outobj.associate_redshift(redshift, redshift_lower, redshift_upper)
    outobj.associate_plot_as(plot_as)
    outobj.associate_cosmology(cosmology)

    output_path = f"{output_directory}/{output_filename}"

    if os.path.exists(output_path):
        os.remove(output_path)

    outobj.write(filename=output_path)

# Also output median and scatter points for composite sample
x_all = np.hstack(x_all) * unyt.dimensionless
y_all = np.hstack(y_all) * unyt.dimensionless

plussig = lambda a: np.percentile(a, 84)
subsig = lambda a: np.percentile(a, 16)

x_bins = np.percentile(x_all, np.linspace(1, 99, 10)) * unyt.dimensionless
x_mids = x_bins[:-1] + np.diff(x_bins) * 0.5

_, y_med, y_sigs = binned_median_line(x_all, y_all, x_bins=x_bins)

y_scatter = unyt.unyt_array((y_sigs[0], y_sigs[1]))
x_scatter = unyt.unyt_array([np.diff(x_bins) * 0.5] * 2)

output_filename = "DeLooze20_composite_median.hdf5"

# Meta-data
comment = f"Median and scatter compiled from multiple surveys"
citation = f"Compiled by De Looze et al. (2020)"
bibcode = "2020MNRAS.496.3668D"
name = "MHI/Mstar as a function of 12 + log10(O/H)"
plot_as = "points"
redshift = 0.0
redshift_lower = 0.0
redshift_upper = 3.0
h = 0.7

# Write everything
outobj = ObservationalData()
outobj.associate_x(
    x_mids, scatter=x_scatter, comoving=True, description="Gas phase 12 + log10(O/H)"
)
outobj.associate_y(
    y_med, scatter=y_scatter, comoving=True, description="HI mass to stellar mass ratio"
)
outobj.associate_citation(citation, bibcode)
outobj.associate_name(name)
outobj.associate_comment(comment)
outobj.associate_redshift(redshift, redshift_lower, redshift_upper)
outobj.associate_plot_as(plot_as)
outobj.associate_cosmology(cosmology)

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

outobj.write(filename=output_path)
