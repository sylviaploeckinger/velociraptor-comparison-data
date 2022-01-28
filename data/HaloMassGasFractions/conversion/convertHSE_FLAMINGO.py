from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys
from astropy.cosmology import FlatLambdaCDM

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Cosmology
h_sim = cosmology.h
Omega_b = cosmology.Ob0
Omega_m = cosmology.Om0
FLATCDM = FlatLambdaCDM(H0=70, Om0=0.30)

output_filename = "HSE-FLAMINGO.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

#The Data
log10_M_500 = [13.6, 13.88571429, 14.05714286, 14.22857143, 14.4, 14.57142857, 14.74285714, 14.91428571]
fb_500 = [0.07378083, 0.08326182, 0.09391293, 0.10494623, 0.11478259, 0.12976777, 0.12958898, 0.13943764]
error_fb_500_p = [0.00371715, 0.00200495, 0.00314589, 0.00490072, 0.00817478, 0.00258178, 0.0021621,  0.00347187]

#Apply Hydrostatic bias
log10_M_500 -=  np.log10(0.74530875)

# Read the data
M_500  = 10**log10_M_500

#Convert to proper units
M_500 = unyt.unyt_array(M_500,units="Msun")
fb_500 = unyt.unyt_array(fb_500,units="dimensionless")
error_fb_500_p = unyt.unyt_array(error_fb_500_p,units="dimensionless")
error_fb_500_m = error_fb_500_p

# Normalise by the cosmic mean
fb_500 = fb_500 / (Omega_b / Omega_m)
error_fb_500_p = error_fb_500_p / (Omega_b / Omega_m)
error_fb_500_m = error_fb_500_m / (Omega_b / Omega_m)

# Define the scatter as offset from the mean value
y_scatter = unyt.unyt_array((error_fb_500_m, error_fb_500_p))

# Meta-data
comment = (
    "Gas fraction data by taking the median for a large set of HSE data. "
    "This is the HSE calibration data for FLAMINGO. Includes fiducial bias. "
    "Data was corrected for the simulation's cosmology."
)
citation = "FLAMINGO HSE Data"
bibcode = " Personal collection"
name = "HSE Gas Fractions"
plot_as = "points"
redshift = 0.1
h = h_sim

# Write everything
processed = ObservationalData()
processed.associate_x(
    M_500, scatter=None, comoving=True, description="Halo mass (M_500)"
)
processed.associate_y(
    fb_500, scatter=y_scatter, comoving=True, description="Gas fraction (<R_500)"
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
