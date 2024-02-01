#!/bin/bash

# read in chips
day=$1
arr=( "$@" )

mkdir data/raw/$day
mkdir configs/$day

# cp configs/template_config.yml configs/$day/
cp -r ~/bucket/data/raw/$day/$day_background/ ./data/raw/$day/

for c in ${arr[@]:1}
  cp -r ~/bucket/data/raw/$day/$c/ ./data/raw/$day/
