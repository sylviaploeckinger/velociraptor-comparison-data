from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

output_directory = "../"
if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# CC. No need to convert Tolstoy+ data to Asplund+, this has already been done when
# preparing the raw data.
# Fe_over_H = 12.0 - 4.5
# Mg_over_H = 12.0 - 4.4
# Mg_over_Fe = Mg_over_H - Fe_over_H

# tabulate/compute the same ratios from Anders & Grevesse (1989)
# Fe_over_H_AG89 = 7.67
# Mg_over_H_AG89 = 7.58

# Mg_over_Fe_AG89 = Mg_over_H_AG89 - Fe_over_H_AG89

for galaxy in ["Carina", "MW", "Fornax", "Sculptor", "Sagittarius"]:
    input_filename = "../raw/{0}.txt".format(galaxy)

    output_filename = "Tolstoy09_{0}.hdf5".format(galaxy)

    data = np.loadtxt(input_filename)
    FeH = data[:, 0]  # + Fe_over_H_AG89 - Fe_over_H
    MgFe = data[:, 1]  # + Mg_over_Fe_AG89 - Mg_over_Fe
    x = unyt.unyt_array(FeH * unyt.dimensionless)
    y = unyt.unyt_array(MgFe * unyt.dimensionless)

    # Meta-data
    comment = (
        "Solar abundances are taken from Asplund et al. (2009), "
        "[Fe/H]Sun = 7.5 and [Mg/H]Sun = 7.6"
    )
    citation = "Tolstoy et al. (2009), {0}".format(galaxy)
    bibcode = "2009ARA&A..47..371T"
    name = "[Mg/Fe] as a function of [Fe/H] for {0}".format(galaxy)
    plot_as = "points"
    redshift = 0.0

    # Write everything
    processed = ObservationalData()
    processed.associate_x(x, scatter=None, comoving=False, description="[Fe/H]")
    processed.associate_y(y, scatter=None, comoving=False, description="[Mg/Fe]")
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
