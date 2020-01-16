cd data

for func in *;
do
    cd $func/conversion

    for convert in ./*.py;
    do
        python3 $convert ../../../cosmology.py;
    done

    cd ../..
done
