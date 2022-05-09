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
    
input_filename = "../raw/dustpedia_cigale_results_final_version.csv"
delimiter = ","

output_directory = "../"
output_filename_werr = "Bianchi2018_Data_werr.hdf5"
output_filename_noerr = "Bianchi2018_Data_noerr.hdf5"


if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data (and convert ids to 0s for now)
converters = {0: lambda s: 0}
raw = np.loadtxt(input_filename, delimiter=delimiter, converters=converters)
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

x_scatter = unyt.unyt_array([M_serr, M_serr])
y_scatter = unyt.unyt_array([M_derr, M_derr])

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star, scatter=x_scatter, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    M_dust, scatter=y_scatter, comoving=True, description="Galaxy Dust Mass"
)
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(redshift, 0, 0.5)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"{output_directory}/{output_filename_werr}"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star, scatter=None, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    M_dust, scatter=None, comoving=True, description="Galaxy Dust Mass"
)
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(redshift, 0, 0.5)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"{output_directory}/{output_filename_noerr}"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)
