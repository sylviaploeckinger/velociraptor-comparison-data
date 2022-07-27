from velociraptor.observations.objects import ObservationalData
from velociraptor.fitting_formulae.smhmr import moster_raw

import unyt
import os
import sys
import numpy as np

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Cosmology
h_sim = cosmology.h

output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)


M_halo, M_star = unyt.unyt_array(
    [
        [56721611.98374207, 12770.33194559625],
        [66697770.34525908, 20268.4996181332],
        [82180533.32338542, 33794.063843929966],
        [97648478.59221399, 55359.3674940112],
        [116533007.44544114, 88143.41033205223],
        [138002489.4983089, 138337.4634628402],
        [165806436.4024939, 216542.39380698395],
        [203461129.85167322, 338639.88757873926],
        [241663678.7559223, 510766.88166467537],
        [301380430.515492, 815815.6218982262],
        [383294215.7194115, 1284168.593466909],
        [539021229.6371884, 2185780.564233959],
        [733929446.0632615, 3614480.4179432406],
        [873778036.1071911, 4204256.880033881],
        [1219588861.6259685, 7335696.761586763],
        [1696231936.852034, 10921579.296876188],
        [2292792334.064888, 17974878.48840379],
        [3240500730.3894515, 26444276.232642855],
        [4257311595.5698156, 37250540.7420971],
        [5714111770.160652, 53288502.198403545],
        [7119576699.547271, 80198907.61316414],
        [8920387100.242449, 128897565.76231846],
        [12306198358.696415, 165320607.50368923],
        [16614314259.340897, 251910368.38282087],
        [21271980617.296864, 395230257.01595145],
        [28145763938.123302, 555176846.9460926],
        [34768505184.97538, 980249882.3529066],
        [44341899574.21948, 1430511343.8168902],
        [58897692434.474686, 1923603624.9569533],
        [73603539415.92117, 2822951378.467494],
        [89491916387.28206, 4273006259.2256036],
        [113630933884.65837, 5215661582.168658],
    ],
    "Solar_Mass",
).T


# Compute the mass ratios and their scatter
MRatio = M_star / M_halo

# Define metadata, which is the same in all the three cases
citation = "Rosdahl et al. (2022)"
bibcode = "2022arXiv220703232R"
name = "The stellar mass-halo mass relation in SPHINX at z=6"
plot_as = "line"
redshift = 6.0

# We purposely make this data show up not only a z=0 but also at higher z
redshift_lower, redshift_upper = 5.0, 7.0
h = h_sim

# stellar mass versus halo mass
comment = (
    "Measurements of stellar mass-halo mass relation from the entire mass "
    "of stars within each halo."
    "Cosmology: Planck 2015 (un-corrected). "
    "The data is extracted by hand from Fig. 4, and is a mean not median. "
)

output_filename = "Rosdahl2022.hdf5"

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_halo,
    scatter=None,
    comoving=False,
    description="Halo Mass ($M_{\\rm vir}$)",
)
processed.associate_y(
    M_star, scatter=None, comoving=False, description="Subhalo Stellar Mass"
)
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(redshift, redshift_lower, redshift_upper)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)

# stellar-to-halo mass ratio versus halo mass
comment = (
    "Measurements of stellar mass-halo mass relation from the entire mass "
    "of stars within each halo."
    "Cosmology: Planck 2015 (un-corrected). "
    "The data is extracted by hand from Fig. 4, and is a mean not median. "
)

output_filename = "Rosdahl2022_Ratio.hdf5"

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_halo,
    scatter=None,
    comoving=False,
    description="Halo Mass ($M_{\\rm vir}$)",
)
processed.associate_y(
    MRatio,
    scatter=None,
    comoving=False,
    description="Subhalo Stellar Mass / Halo Mass ($M_* / M_{\\rm vir}$)",
)
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(redshift, redshift_lower, redshift_upper)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)

# stellar-to-halo mass ratio versus stellar mass
comment = (
    "Measurements of stellar mass-halo mass relation from the entire mass "
    "of stars within each halo."
    "Cosmology: Planck 2015 (un-corrected). "
    "The data is extracted by hand from Fig. 4, and is a mean not median. "
)

output_filename = "Rosdahl2022_RatioStellar.hdf5"

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star, scatter=None, comoving=False, description="Subhalo Stellar Mass"
)
processed.associate_y(
    MRatio,
    scatter=None,
    comoving=False,
    description="Subhalo Stellar Mass / Halo Mass ($M_* / M_{\\rm vir}$)",
)
processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(redshift, redshift_lower, redshift_upper)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)
