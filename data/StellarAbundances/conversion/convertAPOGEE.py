from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys
import h5py

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/APOGEE_data.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

APOGEE_data = h5py.File(input_filename, "r")

element_list = np.array(["C", "MG", "O", "N", "OH", "OHMGFE"])

for element in element_list:

    output_filename = "APOGEE_data_{0}.hdf5".format(element)

    if element == "OH":
        FE_H = APOGEE_data["FE_H"][:]
        O_FE = APOGEE_data["O_FE"][:]
        O_H = O_FE + FE_H
    elif element == "OHMGFE":
        MG_FE = APOGEE_data["MG_FE"][:]
        FE_H = APOGEE_data["FE_H"][:]
        O_FE = APOGEE_data["O_FE"][:]
        O_H = O_FE + FE_H
    else:
        FE_H = APOGEE_data["FE_H"][:]
        el_FE = APOGEE_data[f"{element}_FE"][:]

    # compute COLIBRE assumed abundances ( Asplund et al. 2009 )
    Fe_over_H = 7.5
    Mg_over_H = 7.6
    O_over_H = 8.69
    C_over_H = 8.43
    N_over_H = 7.83

    Mg_over_Fe_AS09 = Mg_over_H - Fe_over_H
    O_over_Fe_AS09 = O_over_H - Fe_over_H
    C_over_Fe_AS09 = C_over_H - Fe_over_H
    N_over_Fe_AS09 = N_over_H - Fe_over_H
    O_over_H_AS09 = O_over_H

    # tabulate/compute the same ratios from Grevesse, Asplund & Sauval (2007)
    Fe_over_H_GA07 = 7.45
    Mg_over_H_GA07 = 7.53
    O_over_H_GA07 = 8.66
    C_over_H_GA07 = 8.39
    N_over_H_GA07 = 7.78

    # --
    Mg_over_Fe_GA07 = Mg_over_H_GA07 - Fe_over_H_GA07
    O_over_Fe_GA07 = O_over_H_GA07 - Fe_over_H_GA07
    C_over_Fe_GA07 = C_over_H_GA07 - Fe_over_H_GA07
    N_over_Fe_GA07 = N_over_H_GA07 - Fe_over_H_GA07

    FE_H += Fe_over_H_GA07 - Fe_over_H
    if element == "O":
        el_FE += O_over_Fe_GA07 - O_over_Fe_AS09
    if element == "MG":
        el_FE += Mg_over_Fe_GA07 - Mg_over_Fe_AS09
    if element == "N":
        el_FE += N_over_Fe_GA07 - N_over_Fe_AS09
    if element == "C":
        el_FE += C_over_Fe_GA07 - C_over_Fe_AS09
    if element == "OH":
        O_H += O_over_H_GA07 - O_over_H_AS09
        O_FE += O_over_Fe_GA07 - O_over_Fe_AS09
    if element == "OHMGFE":
        O_H += O_over_H_GA07 - O_over_H_AS09
        MG_FE += Mg_over_Fe_GA07 - Mg_over_Fe_AS09

    if element == "OH":
        x = O_H
        y = O_FE
    elif element == "OHMGFE":
        x = O_H
        y = MG_FE
    else:
        x = FE_H
        y = el_FE

    x = unyt.unyt_array(x * unyt.dimensionless)
    y = unyt.unyt_array(y * unyt.dimensionless)

    # Meta-data
    comment = "Solar abundances converted to Asplund et al. (2009)"
    citation = "Holtzman, J. A. et al. (2018), {0}".format(element)
    bibcode = "2018AJ....156..125H"
    if element == "OH":
        name = "[O/Fe] as a function of [O/H]".format(element)
        ylabel = "[O/Fe]"
        xlabel = "[O/H]"
    elif element == "OHMGFE":
        name = "[Mg/Fe] as a function of [O/H]".format(element)
        ylabel = "[Mg/Fe]"
        xlabel = "[O/H]"
    else:
        name = "[{0}/Fe] as a function of [Fe/H]".format(element)
        ylabel = "[{0}/Fe]".format(element)
        xlabel = "[Fe/H]"

    plot_as = "points"
    redshift = 0.0

    # Write everything
    processed = ObservationalData()
    processed.associate_x(x, scatter=None, comoving=False, description=xlabel)
    processed.associate_y(y, scatter=None, comoving=False, description=ylabel)
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
