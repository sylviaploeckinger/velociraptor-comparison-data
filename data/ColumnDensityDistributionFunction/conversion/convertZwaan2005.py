from velociraptor.observations.objects import ObservationalData
import unyt
import numpy as np
import os
import sys


def cddf_zwaan():

    # Meta-data
    name = f"CDDF fit from Zwaan et al. (2005)"
    comment = ""

    citation = "Zwaan et al. (2005)"
    bibcode = "2005MNRAS.364.1467Z"
    plot_as = "line"
    output_filename = "Zwaan2005.hdf5"
    output_directory = "../"

    # Create observational data instance
    processed = ObservationalData()
    processed.associate_citation(citation, bibcode)
    processed.associate_name(name)
    processed.associate_comment(comment)
    processed.associate_cosmology(cosmology)

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    # Fit parameters
    NHI_star = 10.0 ** 21.2
    f_star = 0.0193
    beta = 1.24

    logNHI = np.linspace(19.8, 22.1, 100)
    dlogNHI = logNHI[1] - logNHI[0]
    logNHI_minus = logNHI - 0.5 * dlogNHI
    logNHI_plus = logNHI + 0.5 * dlogNHI

    dlogNHI = logNHI_plus - logNHI_minus
    dNHI = 10.0 ** logNHI_plus - 10.0 ** logNHI_minus

    NHI = 10.0 ** logNHI
    f_NHI = (f_star / NHI_star) * (NHI_star / NHI) ** beta * np.exp(-NHI / NHI_star)
    # convert from d/dN to d/dlogN
    f_NHI *= dNHI / dlogNHI

    NHI_bin = unyt.unyt_array(NHI, units="cm**(-2)")
    NHI_scatter = unyt.unyt_array(
        (NHI - 10.0 ** logNHI_minus, 10.0 ** logNHI_plus - NHI),
        units="cm**(-2)",
    )
    f_NHI_bin = unyt.unyt_array(f_NHI, units="dimensionless")

    processed.associate_x(
        NHI_bin, scatter=NHI_scatter, comoving=False, description="Column density"
    )
    processed.associate_y(
        f_NHI_bin,
        scatter=None,
        comoving=False,
        description="Column density distribution function",
    )

    z_minus = 2.0
    z_plus = 5.0
    processed.associate_redshift(0.5 * (z_minus + z_plus), z_minus, z_plus)
    processed.associate_plot_as(plot_as)

    output_path = f"{output_directory}/{output_filename}"

    if os.path.exists(output_path):
        os.remove(output_path)

    processed.write(filename=output_path)


# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Generate the Zwaan et al. (2005) fit
cddf_zwaan()
