#!bash
# Formats all relevant source code files.
# Accepts an extra input argument that is passed to black.
# Useful for --check.

extra_arg=${1}

return_code=0

black cosmology.py $extra_arg
return_code=$(( $return_code > $? ? $return_code : $? ))
black plot_individual_dataset.py $extra_arg
return_code=$(( $return_code > $? ? $return_code : $? ))
black convert.py $extra_arg
return_code=$(( $return_code > $? ? $return_code : $? ))

black data/*/conversion/*.py $extra_arg
return_code=$(( $return_code > $? ? $return_code : $? ))

exit $return_code
