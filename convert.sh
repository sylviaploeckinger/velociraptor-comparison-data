return_code=0

cd data

for func in *;
do
    conversion_path=$func/conversion

    if [ -d $conversion_path ]
    then
        cd $conversion_path

        for convert in ./*.py;
        do
            python3 $convert ../../../cosmology.py;
            this_return_code=$?
            return_code=$(( $return_code > $this_return_code ? $return_code : $this_return_code ))

            if [ $this_return_code -ne 0 ]; then
                echo "Script ${func} -> ${convert} failed."
                echo "Return code: ${this_return_code}."
            fi
        done

        cd ../..
    fi
done

exit $return_code
