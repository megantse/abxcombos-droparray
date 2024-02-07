'''
Description: Make folders and transfer data from VM to bucket
Last updated: 2024 Feb, mtse

'''

import os
import shutil

day = snakemake.input[0]
chips = snakemake.input[1]

os.makedirs('data/raw/'+day, exist_ok=True)
os.makedirs('configs/'+day, exist_ok=True)


shutil.copytree('~/bucket/data/raw/'+day+'/'+day+'$_background/', './data/raw/'+day)

for chip in chips:
    shutil.copytree('~/bucket/data/raw/'+day+'/'+chip+'/', './data/raw/'+day)

shutil.copy('configs/template_config.yml', snakemake.output[0])
