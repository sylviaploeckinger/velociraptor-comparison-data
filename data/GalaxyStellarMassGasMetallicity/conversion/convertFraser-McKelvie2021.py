from velociraptor.observations.objects import ObservationalData

import unyt
import os
import sys
import csv

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Cosmology
h_sim = cosmology.h

input_filename = "../raw/Fraser-McKelvie_2021.csv"
output_filename = "Fraser-McKelvie_2021.hdf5"
output_directory = "../"

Mstar_arr, Zgas_arr = [], []

with open(input_filename, "r") as file:
    data = csv.reader(file, delimiter=",")
    for c, row in enumerate(data):
        if c > 0:

            # The raw stellar mass is given in log10
            Mstar = 10.0 ** float(row[4]),

            # The raw Zgas is 12 + log10 (O/H)
            Zgas = float(row[6])

            Mstar_arr.append(Mstar)
            Zgas_arr.append(Zgas)

Mstar_arr = unyt.unyt_array(Mstar_arr, units="Msun")
Zgas_arr = unyt.unyt_array(Zgas_arr, units="dimensionless")

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Meta-data
comment = (
    "Oxygen-based gas metallicity versus stellar mass for 472 star-forming galaxies "
    "extracted from the SAMI Galaxy Survey. "
    "Data obtained assuming a Chabrier IMF and h=0.704. "
    "The R-calibration of Pilyugin & Grebel (2016) was employed to determine oxygen abundances. "
    "The metallicity is expressed as 12 + log10(O/H). In these units the solar metallicity is 8.69."
)
citation = "Fraser-McKelvie et al. (2021)"
bibcode = " 2021MNRAS.tmp.3132F"
name = "Stellar mass - gas phase metallicity relation "
plot_as = "points"
redshift = 0.1
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    Mstar_arr, scatter=None, comoving=True, description="Galaxy Stellar Mass"
)
processed.associate_y(
    Zgas_arr, scatter=None, comoving=True, description="Gas phase metallicity"
)
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
