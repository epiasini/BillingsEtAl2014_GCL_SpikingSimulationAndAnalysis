# create data folders down to network structure level, and copy the
# network structure files to their locations.

data_dir=/home/ucbtepi/code/network/data
gd_range="1
2
3
4
5
6
7
8
9
10
20"

for gd in $gd_range
do
    origin_dir=$data_dir/network_structures/gd$gd
    target_dir=$data_dir/gmr2.90/gd$gd/s28.74
    mkdir -p $target_dir
    echo copying $origin_dir/GCLconnectivity_full.graphml to $target_dir/gmr2.90_gd"$gd"_s28.74.graphml
    cp $origin_dir/GCLconnectivity_full.graphml $target_dir/gmr2.90_gd"$gd"_s28.74.graphml
done

exit 0
