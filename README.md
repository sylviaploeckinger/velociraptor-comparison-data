Observational Data Collection
=============================

This collection of "observational" (really, comparison - as this
includes some theory datasets as well) are used for comparisons
within the
[VELOCIraptor](https://github.com/swiftsim/velociraptor-python)
library. The repository contains two main things:

1. A selection of plain-text tables extracted from various
   comparison data sets.
2. Accompanying scripts that convert these into h-corrected
   (and otherwise cosmology corrected) values for plotting,
   saving them in HDF5 files.

This functionality is supported by the
[ObservationalData](https://velociraptor-python.readthedocs.io/en/latest/observational_data/index.html)
class within the VELOCIraptor python library. The file format,
and information on how to create a given file, is available
in that documentation.

Contributing
------------

To contribute a dataset to this repository, you will need to
save the data in a plain-text format (including a comment), and
create a conversion script.

This conversion script should be as similar as possible to
the example below:

```python
from velociraptor.observations.objects import ObservationalData

import unyt
import numpy as np
import os
import sys

# Exec the master cosmology file passed as first argument
# These lines are _required_ and you are required to use
# the cosmology specified (this is an astropy.cosmology
# instance)
with open(sys.argv[1], "r") as handle:
    exec(handle.read())

input_filename = "../raw/AuthorYear.txt"
delimiter = None

output_filename = "AuthorYear.hdf5"
output_directory = "../"

if not os.path.exists(output_directory):
    os.mkdir(output_directory)

processed = ObservationalData()
raw = np.loadtxt(input_filename, delimiter=delimiter)

comment = f"Comment About Your Data. h-corrected for SWIFT using Cosmology: {cosmology.name}."
citation = "Author et al. (Year)"
bibcode = "Bibcode from ADS"
name = "Name of Plot"
plot_as = "points/line"
redshift = 0.000000
h = cosmology.h

x = convert_x_to_physical_units
y = convert_y_to_physical_units
# y_scatter should be a 1xN or 2xN array describing offsets from
# the median point 'y'
y_scatter = convert_y_scatter_to_physical_units

processed.associate_x(x, scatter=None, comoving=True, description="x Description")
processed.associate_y(y, scatter=y_scatter, comoving=True, description="y Description")
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
```

There are several examples within the repository already.

You should then organise things as follows:
+ Raw data: `data/PlotName/raw/AuthorYear.txt`
+ Conversion script: `data/PlotName/conversion/convert_AuthorYear.py`

This ensures that it is picked up by the automatic conversion script.

Before creating a pull request, please ensure that you run
the format_all.sh script.
