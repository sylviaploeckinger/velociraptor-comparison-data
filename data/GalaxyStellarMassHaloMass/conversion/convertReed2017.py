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

# Arrays to collect data plots from individual dwarf galaxies
# (Mstar and Mhalo as well as their upper and lower error bars)
M_halo = []
M_halo_err_p = []
M_halo_err_m = []

M_star = []
M_star_err_p = []
M_star_err_m = []

# Conversion factor to solar masses for M_halo (M200)
M_halo_units = 1e10

# Conversion factor to solar masses for M_star
M_star_units = 1e7

# Read the data
with open("../raw/Read2017.txt") as file:
    for line in file.readlines():
        # Skip header
        if line[0] == "#":
            pass
        else:
            # Unpack the line
            _, Ms, Ms_p, Ms_m, _, Mh, Mh_p, Mh_m, _ = line.split(";")

            # Collect M_halo data
            M_halo.append(float(Mh) * M_halo_units)
            M_halo_err_p.append(float(Mh_p) * M_halo_units)
            M_halo_err_m.append(-float(Mh_m) * M_halo_units)

            # Collect M_star data
            M_star.append(float(Ms) * M_star_units)
            M_star_err_p.append(float(Ms_p) * M_star_units)
            M_star_err_m.append(-float(Ms_m) * M_star_units)

# Wrap everything into unyt arrays
M_star = unyt.unyt_array(M_star, units="Msun")
M_star_scatter = unyt.unyt_array((M_star_err_m, M_star_err_p), units="Msun")

M_halo = unyt.unyt_array(M_halo, units="Msun")
M_halo_scatter = unyt.unyt_array((M_halo_err_m, M_halo_err_p), units="Msun")

# Compute the mass ratios and their scatter
MRatio = M_star / M_halo
MRatio_scatter = MRatio * np.sqrt(
    (M_star_scatter / M_star) ** 2 + (M_halo_scatter / M_halo) ** 2
)

# Define metadata, which is the same in all the three cases
citation = "Reed et al. (2017)"
bibcode = "2017MNRAS.467.2019R"
name = "The stellar mass-halo mass relation of isolated field dwarfs."
plot_as = "points"
redshift = 0.0

# We purposely make this data show up not only a z=0 but also at higher z
redshift_lower, redshift_upper = -0.1, 2.1
h = h_sim

comment = (
    "Measurements of stellar mass-halo mass relation obtained by fitting the "
    "rotation curves of isolated dwarf galaxies. "
    "Cosmology: Omega_m=0.27, Omega_lambda=0.73, h=0.7, sigma_8=0.82, "
    "n_s=0.95. (un-corrected). "
    "The data is taken from Table 2 in Reed's paper. "
    "The galaxies whose data was considered as 'bad' (see the discussion on the "
    "exclusion of 'rogues' in paragraph 4.4) are not shown. "
    "The halo mass used is the 'virial' mass within a spherical region whose mean "
    "density is 200 times the critical density of the Universe today. "
    "Plotted is the stellar mass versus halo mass at z=0."
)

output_filename = "Reed2017.hdf5"

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_halo,
    scatter=M_halo_scatter,
    comoving=False,
    description="Halo Mass ($M_{200, {\rm crit}}$)",
)
processed.associate_y(
    M_star,
    scatter=M_star_scatter,
    comoving=False,
    description="Galaxy Stellar Mass",
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

comment = (
    "Measurements of stellar mass-halo mass relation obtained by fitting the "
    "rotation curves of isolated dwarf galaxies. "
    "Cosmology: Omega_m=0.27, Omega_lambda=0.73, h=0.7, sigma_8=0.82, "
    "n_s=0.95. (un-corrected). "
    "The data is taken from Table 2 in Reed's paper. "
    "The galaxies whose data was considered as 'bad' (see the discussion on the "
    "exclusion of 'rogues' in paragraph 4.4) are not shown. "
    "The halo mass used is the 'virial' mass within a spherical region whose mean "
    "density is 200 times the critical density of the Universe today. "
    "Plotted is the stellar-to-halo mass ratio versus halo mass at z=0."
)

output_filename = "Reed2017_Ratio.hdf5"

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_halo,
    scatter=M_halo_scatter,
    comoving=False,
    description="Halo Mass ($M_{200, {\rm crit}}$)",
)
processed.associate_y(
    MRatio,
    scatter=MRatio_scatter,
    comoving=False,
    description="Galaxy Stellar Mass / Halo Mass ($M_{200, {\rm crit}}$)",
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

comment = (
    "Measurements of stellar mass-halo mass relation obtained by fitting the "
    "rotation curves of isolated dwarf galaxies. "
    "Cosmology: Omega_m=0.27, Omega_lambda=0.73, h=0.7, sigma_8=0.82, "
    "n_s=0.95. (un-corrected). "
    "The data is taken from Table 2 in Reed's paper. "
    "The galaxies whose data was considered as 'bad' (see the discussion on the "
    "exclusion of 'rogues' in paragraph 4.4) are not shown. "
    "The halo mass used is the 'virial' mass within a spherical region whose mean "
    "density is 200 times the critical density of the Universe today. "
    "Plotted is the stellar-to-halo mass ratio versus stellar mass at z=0"
)

output_filename = "Reed2017_RatioStellar.hdf5"

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_star,
    scatter=M_star_scatter,
    comoving=False,
    description="Galaxy Stellar Mass",
)
processed.associate_y(
    MRatio,
    scatter=MRatio_scatter,
    comoving=False,
    description="Galaxy Stellar Mass / Halo Mass ($M_{200, {\rm crit}}$)",
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
