from velociraptor.observations.objects import ObservationalData
import re
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
    
input_filename_cigale = "../raw/dustpedia_cigale_results_final_version.csv"
input_filename_h2 = "../raw/dustpedia_H2.csv"
delimiter = ","

output_directory = "../"
output_filename_werr = "DeVis2019_Data_werr.hdf5"
output_filename_noerr = "DeVis2019_Data_noerr.hdf5"


if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data (and convert ids to 0s for now)
converters = {0: lambda s: int.from_bytes(str(s).encode(), 'little')}
raw_cigale = np.loadtxt(input_filename_cigale, delimiter=delimiter, converters=converters)
raw_h2 = np.loadtxt(input_filename_h2, delimiter=delimiter, converters=converters, usecols=range(3))

# x-match catalogues
_, sort1, sort2 = np.intersect1d(raw_cigale[:,0], raw_h2[:,0], return_indices=True)
raw_cigale = raw_cigale[sort1]
raw_h2 = raw_h2[sort2]

M_mol = raw_h2[:, 1] * (unyt.Solar_Mass) * (h_sim / h_obs) ** -2
M_merr = raw_h2[:, 2] * (unyt.Solar_Mass) * (h_sim / h_obs) ** -2
M_dust = raw_cigale[:, 17] * unyt.Solar_Mass * (h_sim / h_obs) ** -2
M_derr = raw_cigale[:, 18] * unyt.Solar_Mass * (h_sim / h_obs) ** -2


lines = []
labels = []
inc = 0

# Meta-data
comment = (
    "Data obtained directly from Dustpedia archive, converting",
    "from Salpeter to Chabrier IMF and using  h=0.7324 (Clark+18). "
    f"h-corrected for SWIFT using cosmology: {cosmology.name}. "
    "M_dust obtained using the CIGALE code (see Bianchi+18), with"
    " HI+H2 values obtained from (De Vis + 2019)"
)
citation = "Dustpedia + CIGALE (Bianchi et. al 2018, De Vis et. al 2019)"
bibcode = "2018A&A...609A..37C"
name = "Molecular Gas Mass - Dust Mass data"
plot_as = "points"
redshift = 0.02
h = h_sim

x_scatter = unyt.unyt_array([M_merr, M_merr])
y_scatter = unyt.unyt_array([M_derr, M_derr])

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_mol, scatter=x_scatter, comoving=True, description="Galaxy Star Formation Rate"
)
processed.associate_y(
    M_dust, scatter=y_scatter, comoving=True, description="Galaxy Molecular Gas Mass"
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
    M_mol, scatter=None, comoving=True, description="Galaxy Molecular Gas Mass"
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
