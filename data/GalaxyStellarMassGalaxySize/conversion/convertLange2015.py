from velociraptor.observations.objects import ObservationalData
from velociraptor.tools.lines import binned_median_line

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = {"H": "../raw/Lange2015HBand.txt", "r": "../raw/Lange2015rBand.txt"}
delimiter = None
half_mass = 1
log_mass = 0

output_filename = {"H": "Lange2015HBand.hdf5", "r": "Lange2015rBand.hdf5"}
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)


for band in ["r", "H"]:
    processed = ObservationalData()
    raw = np.loadtxt(input_filename[band], delimiter=delimiter)

    comment = (
        "Assuming Chabrier IMF (2003), assumes a standard LCDM cosmology with "
        "H=70 km/s/Mpc and O_lambda=0.7. Provided h-free so no h-correction made. "
        f"This data is for the {band}-band only."
    )
    citation = f"Lange et al. (2015, {band}-band)"
    bibcode = "2015MNRAS.447.2603L"
    name = "Galaxy Stellar Mass-Galaxy Size"
    plot_as = "line"
    redshift = 0.05
    h_obs = 0.7
    h = cosmology.h

    half_mass = raw.T[1] * unyt.kpc
    M = 10 ** (raw.T[0]) * unyt.Solar_Mass

    # We now need to bin this data
    mass_bins = unyt.Solar_Mass * np.logspace(7, 12, 26)
    mass_bin_centers, galaxy_sizes, galaxy_size_scatter = binned_median_line(
        x=M, y=half_mass, x_bins=mass_bins
    )

    processed.associate_x(
        mass_bin_centers,
        scatter=None,
        comoving=False,
        description="Galaxy Stellar Mass",
    )
    processed.associate_y(
        galaxy_sizes,
        scatter=galaxy_size_scatter,
        comoving=False,
        description="Galaxy Half-Mass Radius",
    )
    processed.associate_citation(citation, bibcode)
    processed.associate_name(name)
    processed.associate_comment(comment)
    processed.associate_redshift(redshift)
    processed.associate_plot_as(plot_as)
    processed.associate_cosmology(cosmology)

    output_path = f"{output_directory}/{output_filename[band]}"

    if os.path.exists(output_path):
        os.remove(output_path)

    processed.write(filename=output_path)
