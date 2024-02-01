#!/bin/bash

# day number will be $1, copy it over
cd ~/analysis/chips/pilot_screen/data/raw
cp -r ~/bucket/data/$1 .

###
# read in chips
day=$1
arr=( "$@" )

for c in ${arr[@]:1}
  cp -r ~/bucket/data/raw/$1/$c/

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
###
rule vm_transfer:
	input:
		'~/bucket/PRODUCTION/batch1_data/raw/{DAY}/{chip}/', # 0
		'~/bucket/PRODUCTION/batch1_data/raw/{DAY}/{DAY}_background/' # 1
	output:
		'data/raw/{DAY}/{chip}/',
	shell:
		'cp -r {input} {output}'
