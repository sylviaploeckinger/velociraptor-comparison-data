cd data

for func in *;
do
    cd $func

    rm *.hdf5

    cd ..
done
