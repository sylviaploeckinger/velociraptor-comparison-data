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

input_filename1 = f"../raw/data_querejeta_bar.txt"
input_filename2 = f"../raw/data_querejeta_center.txt"
input_filename3 = f"../raw/data_querejeta_disk.txt"
input_filename4 = f"../raw/data_querejeta_interarm.txt"
input_filename5 = f"../raw/data_querejeta_spiral_arms.txt"

processed = ObservationalData()

comment = "Based on the PHANGS sample"
citation = "Querejeta et al. (2021)"
bibcode = "2021A&A...656A.133Q"
name = "Spatially-resolved H2 Gas Surface Density vs Star Formation Rate Surface Density"
plot_as = "points"

# Reading the Ellison 2020 data

array_of_interest = np.arange(-1,3,0.25)
minimum_surface_density = 0.0

sigma_SFR_bar, sigma_H2_bar, sigma_H2_4p3_bar, sigma_H2_N12_bar, sigma_H2_B13_bar = np.genfromtxt(input_filename1, unpack=True, usecols=(7, 9, 10, 11,12), comments="#")

sigma_SFR_center, sigma_H2_center, sigma_H2_4p3_center, sigma_H2_N12_center, sigma_H2_B13_center = np.genfromtxt(input_filename2, unpack=True, usecols=(7, 9, 10, 11,12), comments="#")

sigma_SFR_disk, sigma_H2_disk, sigma_H2_4p3_disk, sigma_H2_N12_disk, sigma_H2_B13_disk = np.genfromtxt(input_filename3, unpack=True, usecols=(7, 9, 10, 11,12), comments="#")

sigma_SFR_interarm, sigma_H2_interarm, sigma_H2_4p3_interarm, sigma_H2_N12_interarm, sigma_H2_B13_interarm = np.genfromtxt(input_filename4, unpack=True, usecols=(7, 9, 10, 11,12), comments="#")

sigma_SFR_spiral_arms, sigma_H2_spiral_arms, sigma_H2_4p3_spiral_arms, sigma_H2_N12_spiral_arms, sigma_H2_B13_spiral_arms = np.genfromtxt(input_filename5, unpack=True, usecols=(7, 9, 10, 11,12), comments="#")

sigma_SFR_total = np.append(sigma_SFR_bar, sigma_SFR_center)
sigma_SFR_total = np.append(sigma_SFR_total, sigma_SFR_disk)
sigma_SFR_total = np.append(sigma_SFR_total, sigma_SFR_interarm)
sigma_SFR_total = np.append(sigma_SFR_total, sigma_SFR_spiral_arms)

sigma_H2_total = np.append(sigma_H2_bar, sigma_H2_center)
sigma_H2_total = np.append(sigma_H2_total, sigma_H2_disk)
sigma_H2_total = np.append(sigma_H2_total, sigma_H2_interarm)
sigma_H2_total = np.append(sigma_H2_total, sigma_H2_spiral_arms)

sigma_SFR_total[~np.isfinite(sigma_SFR_total)] = 1e-20
sigma_SFR_total[sigma_SFR_total<=0] = 1e-20
sigma_H2_total[~np.isfinite(sigma_H2_total)] = 1e-20
sigma_H2_total[sigma_H2_total <=0] = 1e-20

Sigma_H2 = sigma_H2_total/1.36 # a factor of 1.36 to account for heavy elements
Sigma_SFR = sigma_SFR_total 

binned_data = bin_data_general(np.log10(Sigma_H2), np.log10(Sigma_SFR), array_of_interest, minimum_surface_density)

SigmaH2 = unyt.unyt_array(10**binned_data[0], units="Msun/pc**2")

SigmaSFR = unyt.unyt_array(10**binned_data[1], units="Msun/yr/kpc**2")

SigmaSFR_err = unyt.unyt_array([np.abs(10**(binned_data[1]) - 10**(binned_data[1]-binned_data[2])), np.abs(10**(binned_data[1]+binned_data[3]) -10**(binned_data[1]))], units="Msun/yr/kpc**2")

processed.associate_x(SigmaH2, scatter=None, comoving=False, description="H2 Surface density")

processed.associate_y(SigmaSFR, scatter=SigmaSFR_err, comoving=False, description="Star Formation Rate Surface Density")

processed.associate_citation(citation, bibcode)
processed.associate_name(name)
processed.associate_comment(comment)
processed.associate_redshift(0.0, 0.0, 0.0)
processed.associate_plot_as(plot_as)
processed.associate_cosmology(cosmology)

output_path = f"../Querejeta2021.hdf5"

if os.path.exists(output_path):
    os.remove(output_path)

processed.write(filename=output_path)

