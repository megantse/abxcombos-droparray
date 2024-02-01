#!/bin/bash/

# read in chips
day=$1
arr=( "$@" )

# convert to appropriate expected output files:
# set an indicator variable
indicator=1
for c in ${arr[@]:1}
do
	if [[ ! -f "data/output/"$day"/"$c"_qcfiltered.csv" ]]; then
		indicator=0
		echo "file not found: data/output/"$day"/"$c"_qcfiltered.csv"
	fi
done

# if all output files exist, write dummy text file
if ((indicator)); then
	touch "data/output/"$day"/qc_notebook_output.txt"
	echo "pass"
fi
