from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Mishenina_1999.txt"

output_directory = "../"
output_filename = "Mishenina00_data.hdf5"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

Fe_over_H = 12.0 - 4.5
O_over_H = 12.0 - 3.31
O_over_Fe = O_over_H - Fe_over_H

# tabulate/compute the same ratios from Anders & Grevesse (1989)
Fe_over_H_AG89 = 7.67
O_over_H_AG89 = 8.93

O_over_Fe_AG89 = O_over_H_AG89 - Fe_over_H_AG89

data = np.loadtxt(input_filename, skiprows=3, usecols=[1, 3])
FeH_mish = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
OFe_mish = data[:, 1] + O_over_Fe_AG89 - O_over_Fe
x = unyt.unyt_array(FeH_mish * unyt.dimensionless)
y = unyt.unyt_array(OFe_mish * unyt.dimensionless)

# Meta-data
comment = (
    "Solar abundances are taken from Asplund et al. (2009), "
    "[Fe/H]Sun = 7.5 and [Mg/H]Sun = 7.6"
)
citation = "Mishenina et al. (2000), MW"
bibcode = "2000A&A...353..978M"
name = "[O/Fe] as a function of [Fe/H]"
plot_as = "points"
redshift = 0.0

# Write everything
processed = ObservationalData()
processed.associate_x(x, scatter=None, comoving=False, description="[Fe/H]")
processed.associate_y(y, scatter=None, comoving=False, description="[O/Fe]")
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(redshift)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)
