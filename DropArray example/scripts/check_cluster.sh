#!/bin/bash/

# read in chips
day=$1
arr=( "$@" )

# convert to appropriate expected output files:
# set an indicator variable
indicator=1
for c in ${arr[@]:1}
do
	if [[ ! -f "data/interim/"$day"/"$c"_clustered.csv" ]]; then
		indicator=0
		echo "file not found: data/interim/"$day"/"$c"_clustered.csv"
	fi
done

# if all output files exist, write dummy text file
if ((indicator)); then
	touch "data/interim/"$day"/cluster_notebook_output.txt"
	echo "pass"
fi
