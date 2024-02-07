'''
Script: pre_post.py
Last update: 2024 Feb, mtse

Description: mapping pre- and post- merge images
Input:
    1. data/interim/{DAY}/{chip}_clustered.csv
    2. configs/{DAY}/{chip}_config.yml
    3. data/interim/{DAY}/cluster_notebook_output.txt
Output:
    1. data/interim/{DAY}/{chip}_pre_post.csv
    2. data/interim/{DAY}/{chip}_condensed.csv
'''
##############################################################################
################################## IMPORTS ###################################
import pandas as pd
import yaml
import warnings
warnings.filterwarnings("ignore")
import sys
sys.path.insert(1, './')
# kchip imports
import kchip_py3.analyze as kchip_analyze

##############################################################################
################################# PRE POST ###################################
droplets = pd.read_csv(snakemake.input[0],index_col=0)

# read in config file
config_file = snakemake.input[1]
with open(config_file, 'r') as ymlfile:
    config = yaml.safe_load(ymlfile)

# identify premerge wells
pre_wells = droplets.groupby(['IndexX','IndexY','Well_ID'],as_index=False)[['ImageX','ImageY']].mean()

# list of post-merge timepoints
timepoints = set(config['image']['names'].keys())-set(['premerge'])

# analyze data for each timepoint
pre_post_all = []
for timepoint in timepoints:
    # identify postmerge wells and map to pre-merge wells
    pre_post_all.append(kchip_analyze.map_pre_to_post(config,timepoint,pre_wells))

# condense output
condensed = kchip_analyze.stack_timepoints(droplets,pre_post_all,timepoints)

# pre_post output
pre_post = pd.DataFrame(pre_post_all[0])

# save outputs
pre_post.to_csv(snakemake.output[0])
condensed.to_csv(snakemake.output[1])
