# convert data directory structure to new format.
# usage: convert_data_structure.sh

data_dir=/home/ucbtepi/code/network/data

for gd in {1..20}
do
    for mf in 0.{1..9}0
    do
	for cr in {0..1}
	do
	    pattern_dir=$data_dir/gd$gd/cr$cr/iscs0.00/f$mf
	    target_dir_old=$pattern_dir/b0000
	    target_dir_new=$pattern_dir/b1.00
	    mv -v $target_dir_old $target_dir_new
	done
    done
done

exit 0
