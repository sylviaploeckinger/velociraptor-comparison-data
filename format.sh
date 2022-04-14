#!/bin/bash
# Formats all relevant source code files.
# Accepts an extra input argument that is passed to black.
# Useful for --check.

# Check if we can run pip
# This also serves as a check for python3
python3 -m pip --version > /dev/null
if [[ $? -ne 0 ]]
then
  echo "ERROR: cannot run 'python3 -m pip'"
  exit 1
fi

# Check if the virtual environment with black exists
if [ ! -d black_formatting_env ]
then
  echo "Formatting environment not found, installing it..."
  python3 -m venv black_formatting_env
  ./black_formatting_env/bin/python3 -m pip install click==8.0.4 black==21.12b0
fi
# Now we know exactly which black to use
black="./black_formatting_env/bin/python3 -m black"

extra_arg=${1}

return_code=0

$black cosmology.py $extra_arg
return_code=$(( $return_code > $? ? $return_code : $? ))
$black plot_individual_dataset.py $extra_arg
return_code=$(( $return_code > $? ? $return_code : $? ))
$black convert.py $extra_arg
return_code=$(( $return_code > $? ? $return_code : $? ))

$black data/*/conversion/*.py $extra_arg
return_code=$(( $return_code > $? ? $return_code : $? ))

exit $return_code
