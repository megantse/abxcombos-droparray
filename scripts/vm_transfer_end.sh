#!/bin/bash

# read in day number and copy files from VM to bucket
day=$1

cp -r ~/analysis/data/output/$day/*.csv ~/bucket/data/output
cp -r ~/analysis/data/interim/$day/*.csv ~/bucket/data/interim
cp -r ~/analysis/configs/$day/*.yml ~/bucket/configs
cp -r ~/analysis/notebooks/$day/*.ipynb ~/bucket/notebooks
cp -r ~/analysis/figures/$day/*.png ~/bucket/figures
cp -r ~/analysis/figures/$day/*.pdf ~/bucket/figures

rm -r ~/analysis/data/*/$day/
rm -r ~/analysis/*/$day/

sudo shutdown -h now
