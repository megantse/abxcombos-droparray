#!/bin/bash

# read in chips
day=$1

cp -r ~/analysis/data/output/$day/*.csv ~/bucket/data/output
cp -r ~/analysis/data/interim/$day/*.csv ~/bucket/data/interim
cp -r ~/analysis/configs/$day/*.yml ~/bucket/configs
cp -r ~/analysis/notebooks/$day/*.ipynb ~/bucket/notebooks
cp -r ~/analysis/figures/$day/*.png ~/bucket/figures
cp -r ~/analysis/figures/$day/*.pdf ~/bucket/figures

rm -r ~/analysis/data/*/$day/
rm -r ~/analysis/*/$day/

# day number will be $1, copy it over
cd ~/analysis/chips/pilot_screen/data/raw
cp -r ~/bucket/data/$1 .

###
# read in chips
day=$1
arr=( "$@" )

mkdir data/raw/$day
mkdir configs/$day

cp configs/template_config.yml configs/$day/
cp -r ~/bucket/data/raw/$day/$day_background/

for c in ${arr[@]:1}
  cp -r ~/bucket/data/raw/$1/$c/

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


#!/bin/bash

for PLATE in "$@"
do
	cd $PLATE

	mkdir -p input_ph

	cd input_ph

	mkdir -p metadata

	for DATASET in "DAPI-GFP-A594-AF750" "GFP" "A594-AF750" "DAPI-GFP-A594" "DAPI-GFP"
	do
		gsutil -m cp -r "gs://cheeseman-lasagna/${PLATE}/input_ph/${DATASET}" .
	done

	python /luke-perm/scripts/preprocess_plate_ph.py

	gsutil -m cp -r -n preprocess "gs://cheeseman-lasagna/${PLATE}/input_ph"

	gsutil -m cp -r -n metadata "gs://cheeseman-lasagna/${PLATE}/input_ph"

	# for DATASET in "DAPI-GFP-A594-AF750" "GFP" "A594-AF750" "DAPI-GFP-A594" "DAPI-GFP"
	# do
	# 	rm -r ${DATASET}
	# done

	cd ../..

done

sudo shutdown -h now
