#!bash
# Formats all relevant source code files.
# Accepts an extra input argument that is passed to black.
# Useful for --check.

extra_arg=${1}

black cosmology.py $extra_arg
black plot_individual_dataset.py $extra_arg

black data/*/conversion/*.py $extra_arg
