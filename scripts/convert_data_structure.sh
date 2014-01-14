# convert data directory structure to new format.
# usage: convert_data_structure.sh [clean]

data_dir=/home/ucbtepi/code/network/data

for gd in {1..20}
do
    for mf in 0.{1..9}0
    do
	for cr in {0..1}
	do
	    pattern_dir=$data_dir/gd$gd/cr$cr/iscs0.00/f$mf
	    mi_dir_old=$pattern_dir/b0000/sm80/ss0/nm10/ns0
	    mi_dir_new=$pattern_dir/b0000/dta0/mod0/sm80/ss0/nm10/ns0
	    mkdir -p $mi_dir_new
	    mv -v $mi_dir_old/*hdf5 $mi_dir_new
	    case "$1" in
		clean) rm -rv $pattern_dir/b0000/sm80
	    esac
	done
    done
done

exit 0
