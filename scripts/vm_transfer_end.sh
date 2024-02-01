#!/bin/bash

# read in chips
day=$1

cp -r ~/analysis/data/output/$day/*.csv ~/bucket/charlesriver_cs-1/data/output
cp -r ~/analysis/data/interim/$day/*.csv ~/bucket/charlesriver_cs-1/data/interim
cp -r ~/analysis/configs/$day/*.yml ~/bucket/charlesriver_cs-1/configs
cp -r ~/analysis/notebooks/$day/*.ipynb ~/bucket/charlesriver_cs-1/notebooks
cp -r ~/analysis/figures/$day/*.png ~/bucket/charlesriver_cs-1/figures
cp -r ~/analysis/figures/$day/*.pdf ~/bucket/charlesriver_cs-1/figures

rm -r ~/analysis/data/*/$day/
rm -r ~/analysis/*/$day/

sudo shutdown -h now
