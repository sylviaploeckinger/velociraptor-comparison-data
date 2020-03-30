return_code=0

cd data

for func in *;
do
    cd $func/conversion

    for convert in ./*.py;
    do
	echo $convert
        python3 $convert ../../../cosmology.py;
        this_return_code=$?
        return_code=$(( $return_code > $this_return_code ? $return_code : $this_return_code ))

        if [ $this_return_code -ne 0 ]; then
            echo "Script ${func} -> ${convert} failed."
            echo "Return code: ${this_return_code}."
        fi
    done

    cd ../..
done

exit $return_code
