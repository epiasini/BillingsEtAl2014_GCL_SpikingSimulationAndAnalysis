# create data folders down to network structure level, and copy the
# network structure files to their locations.

data_dir=/home/ucbtepi/code/network/data

for gd in {1..20}
do
    origin_dir=$data_dir/network_structures/gd$gd
    mkdir -p $data_dir/gd$gd/cr0
    mkdir -p $data_dir/gd$gd/cr1
    cp -v $origin_dir/GCLconnectivity_full.graphml $data_dir/gd$gd/cr0/gd"$gd"_cr0.graphml
    cp -v $origin_dir/GCLconnectivity_full_randomised.graphml $data_dir/gd$gd/cr1/gd"$gd"_cr1.graphml
done

exit 0
