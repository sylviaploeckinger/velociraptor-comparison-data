from velociraptor.observations.objects import ObservationalData
import unyt
import numpy as np
import os
import sys


# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Meta-data
name = f"Column Density Distribution Function (FIREbox)"
comment = "Curves for all redshift bins [0,6] have been fitted by eye."
citation = "Feldmann et al. (2022, FIREbox)"
bibcode = "2022arXiv220515325F"
plot_as = "line"
output_filename = "FIREbox.hdf5"
output_directory = "../"

# Create observational data instance
processed = ObservationalData()
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_cosmology(cosmology)

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Load raw data
data = np.loadtxt(f"../raw/FIREbox.txt")

logNHI = data[:, 0]
f_NHI = 10.0 ** data[:, 1]

# create NHI bins from the given values
# for each interval in log space, we assign half to the lower bin and
# half to the upper bin. We further assume that the lowest and highest
# bin are symmetric (in log space) around the central value
logNHI_plus = np.zeros(logNHI.shape)
logNHI_minus = np.zeros(logNHI.shape)
logNHI_plus[:-1] = 0.5 * (logNHI[1:] + logNHI[:-1])
logNHI_minus[1:] = 0.5 * (logNHI[1:] + logNHI[:-1])
logNHI_plus[-1] = 2.0 * logNHI[-1] - logNHI_minus[-1]
logNHI_minus[0] = 2.0 * logNHI[0] - logNHI_plus[0]

dlogNHI = logNHI_plus - logNHI_minus
dNHI = 10.0 ** logNHI_plus - 10.0 ** logNHI_minus

# convert from d/dN to d/dlogN
f_NHI *= dNHI / dlogNHI


NHI_bin = unyt.unyt_array(10.0 ** logNHI, units="cm**(-2)")
f_NHI_bin = unyt.unyt_array(f_NHI, units="dimensionless")

processed.associate_x(
    NHI_bin, scatter=None, comoving=False, description="Column density"
)
processed.associate_y(
    f_NHI_bin,
    scatter=None,
    comoving=False,
    description="Column density distribution function",
)

z_minus = 0.0
z_plus = 6.0
processed.associate_redshift(0.5 * (z_minus + z_plus), z_minus, z_plus)
processed.associate_plot_as(plot_as)

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)
