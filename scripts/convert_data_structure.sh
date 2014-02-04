# convert data directory structure to new format.
# usage: convert_data_structure.sh

data_dir=/home/ucbtepi/code/network/data

for gd in {1..20}
do
    for mf in 0.{1..9}0
    do
	for cr in {0..1}
	do
	    for dta in 0.0 0.1 0.3 1.0
	    do
		dta_dir=$data_dir/gd$gd/cr$cr/iscs0.00/f$mf/b1.00/dta$dta
		old_dir=$dta_dir/mod0
		target_dir=$dta_dir/ecs0.0
		mkdir $target_dir
		mv -v $old_dir $target_dir
		done
	done
    done
done

exit 0
