from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/Rodney2014.dat"

processed = ObservationalData()
raw = np.loadtxt(input_filename)

comment = "Based on the CANDELS survey."
citation = "Rodney et al. (2014)"
bibcode = "2014AJ....148...13R"
name = "Cosmic SNIa rate"
plot_as = "points"

z = raw[:, 0]
ratenu = raw[:, 1]
sys_err_p = raw[:, 2]
sys_err_m = raw[:, 3]
stat_p = raw[:, 4]
stat_m = raw[:, 5]

h = cosmology.h

a = 1.0 / (1.0 + z)
err_p = np.sqrt(sys_err_p ** 2 + stat_p ** 2) * (h / 0.7) ** 3
err_m = np.sqrt(sys_err_m ** 2 + stat_m ** 2) * (h / 0.7) ** 3
SNIa_rate = ratenu * 1.0e-5 * (h / 0.7) ** 3
SNIa_err_m = err_m * 1.0e-5
SNIa_err_p = err_p * 1.0e-5

a = unyt.unyt_array(a, units=unyt.dimensionless)
SNIa_rate = unyt.unyt_array(SNIa_rate, units=1.0 / (unyt.yr * unyt.Mpc ** 3))
SNIa_scatter = unyt.unyt_array(
    (SNIa_err_m, SNIa_err_p), units=1.0 / (unyt.yr * unyt.Mpc ** 3)
)

processed.associate_x(
    a, scatter=None, comoving=False, description="Cosmic scale factor"
)
processed.associate_y(
    SNIa_rate,
    scatter=SNIa_scatter,
    comoving=False,
    description="Cosmic SNIa rate",
)
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
zmin = z.min()
zmax = z.max()
processed.associate_redshift(0.5 * (zmin + zmax), zmin, zmax)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"../Rodney2014.hdf5"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)
