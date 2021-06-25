from velociraptor.observations.objects import (
    ObservationalData,
    MultiRedshiftObservationalData,
)
from astropy.cosmology import WMAP9 as cosmology_paper
import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_redshifts = [[0.0, 0.2], [0.5, 0.5], [1, 1], [1.5, 1.5], [2, 2], [3, 3]]
z_list = [0.2, 0.5, 1, 1.5, 2, 3]

output_filename = "Leja_2020.hdf5"
output_directory = "../"
comment = f"Uses panchromatic SED models that infer systematically higher masses and lower star formation rates than standard approaches. h-corrected for SWIFT using cosmology: {cosmology_paper.name}. It suses Chabrier (2003) initial mass function"
citation = "Leja et al. (2020)"
bibcode = "2020ApJ...893..111L"
name = "Continuity GSMF from Leja (2020)"


def schechter(logm, logphi, logmstar, alpha, m_lower=None):
    """
    Generate a Schechter function (in dlogm).
    """
    phi = (
        (10 ** logphi)
        * np.log(10)
        * 10 ** ((logm - logmstar) * (alpha + 1))
        * np.exp(-10 ** (logm - logmstar))
    )
    return phi


def parameter_at_z0(y, z0, z1=0.2, z2=1.6, z3=3.0):
    """
    Compute parameter at redshift ‘z0‘ as a function
    of the polynomial parameters ‘y‘ and the
    redshift anchor points ‘z1‘, ‘z2‘, and ‘z3‘.
    """
    y1, y2, y3 = y
    a = ((y3 - y1) + (y2 - y1) / (z2 - z1) * (z1 - z3)) / (
        z3 ** 2 - z1 ** 2 + (z2 ** 2 - z1 ** 2) / (z2 - z1) * (z1 - z3)
    )
    b = ((y2 - y1) - a * (z2 ** 2 - z1 ** 2)) / (z2 - z1)
    c = y1 - a * z1 ** 2 - b * z1
    return a * z0 ** 2 + b * z0 + c


# Continuity model median parameters + 1-sigma uncertainties.
pars = {
    "logphi1": [-2.44, -3.08, -4.14],
    "logphi1_err": [0.02, 0.03, 0.1],
    "logphi2": [-2.89, -3.29, -3.51],
    "logphi2_err": [0.04, 0.03, 0.03],
    "logmstar": [10.79, 10.88, 10.84],
    "logmstar_err": [0.02, 0.02, 0.04],
    "alpha1": [-0.28],
    "alpha1_err": [0.07],
    "alpha2": [-1.48],
    "alpha2_err": [0.1],
}


def phi_z(z0):
    # Draw samples from posterior assuming independent Gaussian uncertainties.
    # Then convert to mass function at ‘z=z0‘.
    draws = {}
    ndraw = 1000
    for par in ["logphi1", "logphi2", "logmstar", "alpha1", "alpha2"]:
        samp = np.array(
            [
                np.random.normal(median, scale=err, size=ndraw)
                for median, err in zip(pars[par], pars[par + "_err"])
            ]
        )
        if par in ["logphi1", "logphi2", "logmstar"]:
            draws[par] = parameter_at_z0(samp, z0)
        else:
            draws[par] = samp.squeeze()

    # Generate Schechter functions.
    logm = np.linspace(8, 12, 100)[:, None]  # log(M) grid
    phi1 = schechter(
        logm, draws["logphi1"], draws["logmstar"], draws["alpha1"]  # primary component
    )
    phi2 = schechter(
        logm,
        draws["logphi2"],  # secondary component
        draws["logmstar"],
        draws["alpha2"],
    )
    phi = phi1 + phi2  # combined mass function

    # Compute median and 1-sigma uncertainties as a function of mass.
    phi_50, phi_84, phi_16 = np.percentile(phi, [50, 84, 16], axis=1)

    return logm.flatten(), phi_50, phi_84, phi_16


if not os.path.exists(output_directory):
    os.mkdir(output_directory)

multi_z = MultiRedshiftObservationalData()
multi_z.associate_citation(citation, bibcode)
multi_z.associate_name(name)
multi_z.associate_comment(comment)
multi_z.associate_cosmology(cosmology)
multi_z.associate_maximum_number_of_returns(1)

for z, redshifts in zip(z_list, input_redshifts):
    processed = ObservationalData()

    plot_as = "line"
    redshift = z
    redshift_lower, redshift_upper = redshifts
    h_paper = 0.697
    h = cosmology.h

    m, phi, err_p, err_m = phi_z(z)

    log_M = m + 2 * np.log10(h_paper / h)
    M = 10 ** (log_M) * unyt.Solar_Mass
    Phi = phi * (h / h_paper) ** 3 * unyt.Mpc ** (-3)

    processed.associate_x(
        M, scatter=None, comoving=True, description="Galaxy Stellar Mass"
    )
    processed.associate_y(Phi, scatter=None, comoving=True, description="Phi (GSMF)")
    processed.associate_redshift(redshift, redshift_lower, redshift_upper)
    processed.associate_plot_as(plot_as)

    multi_z.associate_dataset(processed)

output_path = f"{output_directory}/{output_filename}"

if os.path.exists(output_path):
    os.remove(output_path)

multi_z.write(filename=output_path)
