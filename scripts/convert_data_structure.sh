# convert data directory structure to new format.

data_dir=/home/ucbtepi/code/network/data

for gd in {1..20}
do
    for mf in 0.{1..9}0
    do
	pattern_dir_old=$data_dir/gmr2.90/gd$gd/s28.74/f$mf
	pattern_dir_new=$data_dir/gd$gd/cr0/iscs0.00/f$mf
	mkdir -p $pattern_dir_new
	cp -v $pattern_dir_old/gmr2.90_gd"$gd"_s28.74_sp128_stim.txt  $pattern_dir_new/sp128_stim.txt
	mi_dir_old=$pattern_dir_old/b00/sm80/ss0/nm10/ns0
	mi_dir_new=$pattern_dir_new/b0000/sm80/ss0/nm10/ns0
	mkdir -p $mi_dir_new
	cp -v $mi_dir_old/sp128_t50.hdf5 $mi_dir_new/sp128_t50_sdur150.hdf5
    done
done

exit 0
