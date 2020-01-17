return_code=0

cd data

for func in *;
do
    cd $func/conversion

    for convert in ./*.py;
    do
        python3 $convert ../../../cosmology.py;
        this_return_code=$?
        return_code=$(( $return_code > $this_return_code ? $return_code : $this_return_code ))
    done

    cd ../..
done

exit $return_code
