from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Trujillo2020_full.txt"
delimiter = None
half_mass = 11
R_1 = 14
log_mass = 15

half_mass_output_filename = "Trujillo2020HalfMass.hdf5"
R1_mass_output_filename = "Trujillo2020R1.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

processed = ObservationalData()
raw = np.loadtxt(
    input_filename, delimiter=delimiter, usecols=[half_mass, R_1, log_mass]
)

comment = f"Assuming Chabrier IMF (2003), and a standard LCDM cosmology with H=70km/s/Mpc and O_L = 0.7. h-corrected for SWIFT using cosmology: {cosmology.name}."
citation = "Trujillo et al. (2020)"
bibcode = "2020arXiv200102689T"
name = "Galaxy Stellar Mass-Galaxy Size"
plot_as = "points"
redshift = 0.0
h_obs = 0.7
h = cosmology.h

half_mass = raw.T[0] * unyt.kpc
R_1 = raw.T[1] * unyt.kpc
M = 10 ** (raw.T[2]) * unyt.Solar_Mass

for radius, output_filename, desc in zip(
    [R_1, half_mass],
    [R1_mass_output_filename, half_mass_output_filename],
    [
        "R_1, radius at which stellar density is 1 msun / pc^2",
        "Galaxy half-mass radius",
    ],
):
    processed.associate_x(
        M, scatter=None, comoving=False, description="Galaxy Stellar Mass"
    )
    processed.associate_y(radius, scatter=None, comoving=False, description=desc)
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
