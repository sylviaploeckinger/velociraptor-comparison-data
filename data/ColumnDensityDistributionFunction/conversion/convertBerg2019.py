from velociraptor.observations.objects import ObservationalData
import unyt
import numpy as np
import os
import sys


def cddf_berg():

    # Meta-data
    comment = ""

    citation = "Berg et al. (2019)"
    bibcode = "2019MNRAS.488.4356B"
    plot_as = "points"
    output_directory = "../"

    for dataset in ["subDLA", "DLA"]:

        output_filename = f"Berg2019_{dataset}.hdf5"
        name = f"CDDF from Berg et al. (2019) - {dataset}s"

        # Create observational data instance
        processed = ObservationalData()
        processed.associate_citation(citation, bibcode)
        processed.associate_name(name)
        processed.associate_comment(comment)
        processed.associate_cosmology(cosmology)

        if not os.path.exists(output_directory):
            os.mkdir(output_directory)

        # Load raw data
        data = np.loadtxt(f"../raw/cddf_Berg2019_{dataset}.dat")

        # Fetch the fields we need
        logNHI = data[:, 0]
        f_NHI = data[:, 3]
        f_NHI_minus = data[:, 2]
        f_NHI_plus = data[:, 4]

        # create NHI bins from the given values
        # for each interval in log space, we assign half to the lower bin and
        # half to the upper bin. We further assume that the lowest and highest
        # bin are symmetric (in log space) around the central value
        logNHI_plus = np.zeros(logNHI.shape)
        logNHI_minus = np.zeros(logNHI.shape)
        logNHI_plus[:-1] = 0.5 * (logNHI[1:] + logNHI[:-1])
        logNHI_minus[1:] = 0.5 * (logNHI[1:] + logNHI[:-1])
        logNHI_plus[-1] = 2.0 * logNHI[-1] - logNHI_minus[-1]
        logNHI_minus[0] = 2.0 * logNHI[0] - logNHI_plus[0]

        dlogNHI = logNHI_plus - logNHI_minus
        dNHI = 10.0 ** logNHI_plus - 10.0 ** logNHI_minus

        # convert from d/dN to d/dlogN
        f_NHI *= dNHI / dlogNHI
        f_NHI_minus *= dNHI / dlogNHI
        f_NHI_plus *= dNHI / dlogNHI

        NHI_bin = unyt.unyt_array(10.0 ** logNHI, units="cm**(-2)")
        NHI_scatter = unyt.unyt_array(
            (
                10.0 ** logNHI - 10.0 ** logNHI_minus,
                10.0 ** logNHI_plus - 10.0 ** logNHI,
            ),
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

        z_minus = 2.3
        z_plus = 3.2
        processed.associate_redshift(0.5 * (z_minus + z_plus), z_minus, z_plus)
        processed.associate_plot_as(plot_as)

        output_path = f"{output_directory}/{output_filename}"

        if os.path.exists(output_path):
            os.remove(output_path)

        processed.write(filename=output_path)


# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Generate, format and save the Berg et al. (2019) data
cddf_berg()
