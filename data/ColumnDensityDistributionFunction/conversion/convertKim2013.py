from velociraptor.observations.objects import ObservationalData
import unyt
import numpy as np
import os
import sys


def cddf_kim():

    # Meta-data
    name = f"CDDF from Kim et al. (2013)"
    comment = ""

    citation = "Kim et al. (2013)"
    bibcode = "2013A&A...552A..77K"
    plot_as = "points"
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
    data = np.loadtxt(f"../raw/cddf_Kim2013.dat")

    # Fetch the fields we need
    logNHI = data[:, 0]

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

    for zm, zp, ofs in zip([1.9, 1.9, 2.4], [3.2, 2.4, 3.2], [1, 4, 7]):

        output_filename = f"Kim2013_z{zm:.1f}_{zp:.1f}.hdf5"

        f_NHI = data[:, ofs]
        # mask out empty rows
        mask = f_NHI < 0
        f_NHI_plus = data[:, ofs + 1]
        f_NHI_minus = data[:, ofs + 2]
        # if Delta(f)_minus is not gives, assume the same as Delta(f)_plus
        f_NHI_minus[f_NHI_minus == 0] = f_NHI_plus[f_NHI_minus == 0]

        f_NHI_plus = 10.0 ** (f_NHI + f_NHI_plus)
        f_NHI_minus = 10.0 ** (f_NHI - f_NHI_minus)
        f_NHI = 10.0 ** f_NHI

        # convert from d/dN to d/dlogN
        f_NHI *= dNHI / dlogNHI
        f_NHI_minus *= dNHI / dlogNHI
        f_NHI_plus *= dNHI / dlogNHI

        NHI_bin = unyt.unyt_array(10.0 ** logNHI[mask], units="cm**(-2)")
        NHI_scatter = unyt.unyt_array(
            (
                10.0 ** logNHI[mask] - 10.0 ** logNHI_minus[mask],
                10.0 ** logNHI_plus[mask] - 10.0 ** logNHI[mask],
            ),
            units="cm**(-2)",
        )
        f_NHI_bin = unyt.unyt_array(f_NHI[mask], units="dimensionless")
        f_NHI_scatter = unyt.unyt_array(
            (f_NHI[mask] - f_NHI_minus[mask], f_NHI_plus[mask] - f_NHI[mask]),
            units="dimensionless",
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

        z_minus = zm
        z_plus = zp
        processed.associate_redshift(0.5 * (z_minus + z_plus), z_minus, z_plus)
        processed.associate_plot_as(plot_as)

        output_path = f"{output_directory}/{output_filename}"

        if os.path.exists(output_path):
            os.remove(output_path)

        processed.write(filename=output_path)


# Exec the master cosmology file passed as first argument
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

# Generate, format and save the Kim et al. (2013) data
cddf_kim()
