from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys


def stringify_z(z):
    """
    Eagle-style text formatting of redshift label.
    Example: z=1.5 will be printed as z001p500.

    z: The redshift to produce a label for
    """
    whole = int(z)
    frac = int(1000 * (z - whole))
    return f"z{whole:03d}p{frac:03d}"


# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/DaCunha2015.txt"
delimiter = None

output_filename = "DaCunha2015_{}.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

processed = ObservationalData()
raw = np.loadtxt(input_filename, delimiter=delimiter, usecols=range(1, 10))

comment = (
    "Assuming Chabrier IMF. z=0.01 - 0.2. No h-correction "
    "was required as data was supplied h-free. "
)
citation = "Da Cunha et al. (2015) (SMGs, ALESS)"
bibcode = "10.1088/0004-637X/806/1/110"
name = "Galaxy Stellar Mass - Galaxy Dust Mass from ALESS"
plot_as = "points"
h = 0.7

z = raw.T[0]
z_hi = raw.T[0] + raw.T[1]
z_lo = raw.T[0] + raw.T[2]

M = unyt.unyt_array(10 ** raw.T[3], units=unyt.Solar_Mass)
M_hi = unyt.unyt_array(10 ** (raw.T[3] + raw.T[4]), units=unyt.Solar_Mass) - M
M_lo = M - unyt.unyt_array(10 ** (raw.T[3] + raw.T[5]), units=unyt.Solar_Mass)

Mdust = unyt.unyt_array(10 ** raw.T[6], units=unyt.Solar_Mass)
Mdust_hi = unyt.unyt_array(10 ** (raw.T[6] + raw.T[7]), units=unyt.Solar_Mass) - Mdust
Mdust_lo = Mdust - unyt.unyt_array(10 ** (raw.T[6] + raw.T[8]), units=unyt.Solar_Mass)

zbins = np.arange(1, 7) + 0.5
zcens = zbins[:-1] + 0.5 * np.diff(zbins)

# for redshift in zbins:
#     z in

zindex = np.digitize(z, zbins[1:], right=True)


for i in range(zcens.size):
    bdx = zindex == i
    output_path = f"{output_directory}/{output_filename.format(stringify_z(zcens[i]))}"
    print(output_path)

    processed.associate_x(
        M[bdx],
        scatter=unyt.unyt_array((M_lo[bdx], M_hi[bdx])),
        comoving=False,
        description="Galaxy Stellar Mass",
    )
    processed.associate_y(
        Mdust[bdx],
        scatter=unyt.unyt_array((Mdust_lo[bdx], Mdust_hi[bdx])),
        comoving=False,
        description="Galaxy Dust Mass",
    )
    processed.associate_citation(citation, bibcode)
    processed.associate_name(name)
    processed.associate_comment(comment)
    processed.associate_redshift(zcens[i])
    processed.associate_plot_as(plot_as)
    processed.associate_cosmology(cosmology)

    if os.path.exists(output_path):
        os.remove(output_path)

    processed.write(filename=output_path)
