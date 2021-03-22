from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Cosmologies
h_sim = cosmology.h

input_filename = "../raw/Hunter2021_Oh2015_composite.txt"
delimiter = "\t"

output_filename = "Hunter2021_Oh2015_Data.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Read the data
raw = np.loadtxt(input_filename, delimiter=delimiter)
M_star_kin = 1e7 * raw[:, 0] * unyt.Solar_Mass
M_star_sed = 1e7 * raw[:, 1] * unyt.Solar_Mass

R_D = raw[:, 2]
SFR = (pow(10.0, raw[:, 3]) * np.pi * pow(R_D, 2.0)) * unyt.Solar_Mass / unyt.yr

# index galaxies with defined Mstar
M_star = M_star_kin
M_star[M_star_kin < -99.9] = M_star_sed[M_star_kin < -99]
index_defined = np.logical_and(M_star > -99, raw[:, 3] > -99)

M_star = M_star[index_defined]
SFR = SFR[index_defined]

# Meta-data
comment = (
    "Stellar Masses obtained kinematically where available"
    "Otherwise SED stellar masses are used. Measurements  "
    "are cosmology independent as decoupled from the Hubble"
    " flow so remove influence of h. IMF unspecified, though"
    " only a couple of measurements are IMF depenfent (SED "
    "fits). SFRs inferred from FUV luminosities assuming "
    "Chabrier IMF"
)

citation = "LITTLE THINGS (from SFR, Oh et al. 2015, Hunter et al. 2021)"
bibcode = "2021AJ....161...71H"
name = "Stellar mass - Star Formation Rate"
plot_as = "points"
redshift = 0.0
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star, scatter=None, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    SFR, scatter=None, comoving=True, description="Star Formation Rate"
)
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(redshift, 0, 2)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)
