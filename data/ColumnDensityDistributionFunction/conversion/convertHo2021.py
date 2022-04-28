from velociraptor.observations.objects import ObservationalData
import unyt
import numpy as np
import os
import sys


def cddf_ho():

    # Meta-data
    name = f"CDDF from Ho et al. (2021)"
    comment = ""

    citation = "Ho et al. (2021)"
    bibcode = "2021MNRAS.507..704H"
    plot_as = "points"
    output_filename = "Ho2021.hdf5"
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
    data = np.loadtxt(f"../raw/cddf_Ho2021.dat")

    # Fetch the fields we need
    logNHI_minus = data[:, 0]
    logNHI_plus = data[:, 1]
    # add pre-factor 10^{-21}
    f_NHI = data[:, 2] * 1.0e-21
    f_NHI_minus = data[:, 3] * 1.0e-21
    f_NHI_plus = data[:, 4] * 1.0e-21

    logNHI = 0.5 * (logNHI_minus + logNHI_plus)
    dlogNHI = logNHI_plus - logNHI_minus
    dNHI = 10.0 ** logNHI_plus - 10.0 ** logNHI_minus

    # convert from d/dN to d/dlogN
    f_NHI *= dNHI / dlogNHI
    f_NHI_minus *= dNHI / dlogNHI
    f_NHI_plus *= dNHI / dlogNHI

    NHI_bin = unyt.unyt_array(10.0 ** logNHI, units="cm**(-2)")
    NHI_scatter = unyt.unyt_array(
        (10.0 ** logNHI - 10.0 ** logNHI_minus, 10.0 ** logNHI_plus - 10.0 ** logNHI),
        units="cm**(-2)",
    )
    f_NHI_bin = unyt.unyt_array(f_NHI, units="dimensionless")
    f_NHI_scatter = unyt.unyt_array(
        (f_NHI - f_NHI_minus, f_NHI_plus - f_NHI), units="dimensionless"
    )

    processed.associate_x(
        NHI_bin, scatter=NHI_scatter, comoving=False, description="Column density"
    )
    processed.associate_y(
        f_NHI_bin,
        scatter=f_NHI_scatter,
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

# Generate, format and save the Ho et al. (2021) data
cddf_ho()
