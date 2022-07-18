from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

def bin_data_general(array_x, array_y, array_x_bin, x_limit, print_stuff=False):
    # create a SFR value array
    y_array_bin = np.zeros(len(array_x_bin)-1)
    y_array_bin_std_up = np.zeros(len(array_x_bin)-1)
    y_array_bin_std_down = np.zeros(len(array_x_bin)-1)

    for i in range(0,len(array_x_bin)-1):
        mask = (array_x > array_x_bin[i]) & (array_x < array_x_bin[i+1])
        y_array_bin[i] = np.nanmedian(array_y[mask])
        if print_stuff:
            print(array_x_bin[i])
            print(array_y[mask])
            print(stats.describe(array_y[mask]))
        try:
            y_array_bin_std_up[i], y_array_bin_std_down[i] = np.transpose(np.percentile(array_y[mask], [16,84]))
        except:
            y_array_bin_std_up[i], y_array_bin_std_down[i] = [0., 0.]

    array_x_bin =  (array_x_bin[1:] + array_x_bin[:-1])/2.
    y_array_bin_std_up = np.abs(y_array_bin_std_up - y_array_bin)
    y_array_bin_std_down = np.abs(y_array_bin_std_down - y_array_bin)
    mask = array_x_bin > x_limit

    return array_x_bin[mask], y_array_bin[mask], y_array_bin_std_up[mask], y_array_bin_std_down[mask]


# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = f"../raw/Bigiel2008.txt"

processed = ObservationalData()

comment = "Inner regions of galaxies"
citation = "Bigiel et al. (2008)"
bibcode = "2008AJ....136.2846B"
name = "Spatially-resolved H2 + HI Gas Surface Density vs Star Formation Rate Surface Density"
plot_as = "points"

# Reading the bigiel 2008 data

with open(input_filename) as f:
    lines = f.readlines()

size_array = len(lines) - 49

sigma_HI = -5*np.ones(size_array)
sigma_HI_err = -5*np.ones(size_array)
sigma_H2 = -5*np.ones(size_array)
sigma_H2_err = -5*np.ones(size_array)
sigma_SFR = -5*np.ones(size_array)
sigma_SFR_err = -5*np.ones(size_array)

for i in range(49, len(lines)):
    k = i - 49
    word1 = lines[i][25:29]
    word2 = lines[i][30:34]
    word3 = lines[i][35:39]
    word4 = lines[i][40:44]
    word5 = lines[i][45:50]
    word6 = lines[i][52:56]
    if word1 != "    ":
        sigma_HI[k] = float(word1)
    if word2 != "    ":
        sigma_HI_err[k] = float(word2)
    if word3 != "    ":
        sigma_H2[k] = float(word3)
    if word4 != "    ":
        sigma_H2_err[k] = float(word4)
    if word5 != "     ":
        sigma_SFR[k] = float(word5)
    if word6 != "    ":
        sigma_SFR_err[k] = float(word6)

sigma_gas = (10**sigma_H2 + 10**sigma_HI)/1.36 # a factor of 1.36 to account for heavy elements

array_of_interest = np.arange(-1,3,0.25)
minimum_surface_density = 0.4

Obs_Hneutral = sigma_gas 

Obs_SFR = 10**sigma_SFR

binned_data = bin_data_general(np.log10(Obs_Hneutral), np.log10(Obs_SFR), array_of_interest, minimum_surface_density)

SigmaHneutral = unyt.unyt_array(10**binned_data[0], units="Msun/pc**2")

SigmaSFR = unyt.unyt_array(10**binned_data[1], units="Msun/yr/kpc**2")

SigmaSFR_err = unyt.unyt_array([np.abs(10**(binned_data[1]) - 10**(binned_data[1]-binned_data[2])), np.abs(10**(binned_data[1]+binned_data[3]) -10**(binned_data[1]))], units="Msun/yr/kpc**2")

processed.associate_x(SigmaHneutral, scatter=None, comoving=False, description="H2 + HI Surface density")

processed.associate_y(SigmaSFR, scatter=SigmaSFR_err, comoving=False, description="Star Formation Rate Surface Density")

processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(0.0, 0.0, 0.0)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"../Bigiel2008.hdf5"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)

