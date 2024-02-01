import os
import shutil
import getopt, sys

# day = snakemake.input[0]
# chips = snakemake.input[1]

# Get full command-line arguments
full_cmd_arguments = sys.argv

# Keep all but the first
day = full_cmd_arguments[1]
chips = full_cmd_arguments[2]

os.makedirs('data/raw/'+day, exist_ok=True)
os.makedirs('configs/'+day, exist_ok=True)


shutil.copytree('~/bucket/data/raw/'+day+'/'+day+'$_background/', './data/raw/'+day)

for chip in chips:
    shutil.copytree('~/bucket/data/raw/'+day+'/'+chip+'/', './data/raw/'+day)

shutil.copy('configs/template_config.yml', snakemake.output[0])
