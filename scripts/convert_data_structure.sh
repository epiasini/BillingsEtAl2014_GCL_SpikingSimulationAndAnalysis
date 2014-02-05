# convert data directory structure to new format.
# usage: convert_data_structure.sh

data_dir=/home/ucbtepi/code/network/data

for gd in {1..20}
do
    for mf in 0.{1..9}0
    do
	for cr in {0..1}
	do
	    for ecs in 0.0 1.0
	    do
		for dta in 0.0 0.1 0.3 1.0
		do
		    if [ $dta == 0.0 ]; then
			ics=0.0
		    else
			ics=1.0
		    fi
		    dta_dir=$data_dir/gd$gd/cr$cr/iscs0.00/f$mf/b1.00/dta$dta
		    old_dir=$dta_dir/ecs$ecs
		    target_dir=$dta_dir/ics$ics
		    if [ -d "$old_dir" ]; then
			mkdir $target_dir
			mv -v $old_dir $target_dir
		    fi
		done
	    done
	done
    done
done

exit 0
