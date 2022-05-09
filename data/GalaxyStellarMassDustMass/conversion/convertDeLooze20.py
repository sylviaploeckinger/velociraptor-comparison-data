from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys
import copy
import re
from velociraptor.tools.lines import binned_median_line

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Cosmologies
h_obs = 0.7324
h_sim = cosmology.h
imfcorr = 0.58 
    
input_filename_dust = "../raw/dustpedia_cigale_results_final_version.csv"
delimiter = ","

output_directory = "../"
output_filename = "Bianchi2018_Data.hdf5"


if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data
raw = np.loadtxt(input_filename, delimiter=delimiter)
M_star = raw[:, 3] * unyt.Solar_Mass * imfcorr * (h_sim / h_obs) ** -2
M_serr = raw[:, 4] * unyt.Solar_Mass * imfcorr *(h_sim / h_obs) ** -2
M_dust = raw[:, 17] * unyt.Solar_Mass * (h_sim / h_obs) ** -2
M_derr = raw[:, 18] * unyt.Solar_Mass * (h_sim / h_obs) ** -2


lines = []
labels = []
inc = 0

# Meta-data
comment = (
    "Data obtained directly from Dustpedia archive, converting",
    "from Salpeter to Chabrier IMF and using  h=0.7324 (Clark+18). "
    f"h-corrected for SWIFT using cosmology: {cosmology.name}. "
    "Values obtained using the CIGALE code (see Bianchi+18)"
)
citation = "Dustpedia + CIGALE (Bianchi et. al 2018)"
bibcode = "2018A&A...609A..37C"
name = "Stellar mass - Dust Mass data"
plot_as = "points"
redshift = 0.02
h = h_sim

x_scatter = uny.unyt_array([M_serr, M_serr])
y_scatter = uny.unyt_array([M_derr, M_derr])

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star, scatter=x_scatter, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    Z_star, scatter=y_scatter, comoving=True, description="Galaxy Dust Mass"
)
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(redshift, 0, 0.5)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)


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
    x_all.append(pow(10, raw[:, 1]))

    logMHIMstar_med = pow(10, raw[:, 0]) * unyt.dimensionless
    logMHIMstar_lo = pow(10, raw[:, 2]) * unyt.dimensionless
    logMHIMstar_hi = pow(10, raw[:, 4]) * unyt.dimensionless
    logD2Z_med = pow(10, raw[:, 1]) * unyt.dimensionless
    logD2Z_lo = pow(10, raw[:, 3]) * unyt.dimensionless
    logD2Z_hi = pow(10, raw[:, 5]) * unyt.dimensionless

    # Define the scatter as offset from the mean value
    y_scatter = unyt.unyt_array(
        (logMHIMstar_med - logMHIMstar_lo, logMHIMstar_hi - logMHIMstar_med)
    )
    x_scatter = unyt.unyt_array((logD2Z_med - logD2Z_lo, logD2Z_hi - logD2Z_med))

    # Meta-data
    comment = f"values for individual galaxies derived from {label} survey"
    citation = f"{label} compiled by De Looze et al. (2020)"
    bibcode = "2020MNRAS.496.3668D"
    name = "MHI/Mstar as a function of dust-to-metal ratio"
    plot_as = "points"
    redshift = 0.0
    redshift_lower = 0.0
    redshift_upper = 3.0
    h = 0.7

    # Write everything
    outobj = ObservationalData()
    outobj.associate_x(
        logD2Z_med,
        scatter=x_scatter,
        comoving=True,
        description="Gas phase dust-to-metal ratio",
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
name = "MHI/Mstar as a function of dust-to-metal ratio"
plot_as = "points"
redshift = 0.0
redshift_lower = 0.0
redshift_upper = 3.0
h = 0.7

# Write everything
outobj = ObservationalData()
outobj.associate_x(
    x_mids,
    scatter=x_scatter,
    comoving=True,
    description="Gas phase dust-to-metal ratio",
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
