'''
Script: find_drops.py
Last update: 2021 March, mzhu

Description: DETECTS DROPLETS
Input:
    1. configs/{DAY}/{chip}_config.yml
    2. data/raw/{DAY}/{DAY}_background/background_image.tif ** INCORRECT **
Output:
    1. data/interim/{DAY}/{chip}_droplets_a.csv # droplets found (unprocessed)
    2. data/interim/{DAY}/{chip}_droplets_b.csv # overlap removed, assigned to wells / hashed
    3. data/interim/{DAY}/{chip}_droplets_c.csv # at least 2 droplets per well
    4. data/interim/{DAY}/{chip}_droplets_d.csv # spectral overlap corrected
'''
##############################################################################
################################## IMPORTS ###################################
import yaml
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')
from scipy.spatial.distance import pdist, squareform
import sys
sys.path.insert(1, './')
# kchip imports
import kchip_py3.io as kchip_io
import kchip_py3.analyze as kchip_analyze

def droplets_in_same_well(positions):
    '''Identify droplets that are in the same well based on distance threshold
    Inputs:
    - Positions, (k x 2 ndarray), [x, y] positions of droplet centers
    - close_threshold, max distance between droplet centers to be considred in same well (default = 25)
    Outputs:
    - list of tuples (index 1, index 2) of indexes in positions of droplets in the same well
    '''

    pair_distances = squareform(pdist(positions))
    pair_distances[pair_distances==0]=float('inf')
    ## OLD METHOD
    # close_threshold = int(160*config['image']['objective']/config['image']['bins']/config['image']['pixel_size']+10)
    # droplet1,droplet2 = np.where(pair_distances<close_threshold)
    # zip_list = zip(droplet1,droplet2)

    droplet1 = np.asarray(range(pair_distances.shape[0]))
    droplet2 = np.argmin(pair_distances,axis=0)
    zip_list = zip(droplet1,droplet2)

#     Sort list so small index is always first in tuple, and remove duplicates
    return list(set([tuple(sorted(tup_)) for tup_ in zip_list]))

def remove_overlap(config,df,show=0):
    '''
    Remove droplets found in regions of image that overlap

    Inputs:
        :config: config dictionary
        :df: pandas dataFrame

    Outputs: copy of dataframe, with dropped rows
    '''
    df = df.copy()
    maxX = df['IndexX'].max()
    maxY = df['IndexY'].max()

    overlap = config['image']['size']*(1-config['image']['overlap'])

    rmv_index = ((df['ImageX'] > overlap) & ~(df['IndexX']==maxX)) | ((df['ImageY'] > overlap) & ~(df['IndexY']==maxY))

    return df.drop(df.index[rmv_index]).reset_index(drop=True)

if __name__ == '__main__':

    config_file = snakemake.input[0]
    with open(config_file, 'r') as ymlfile:
        config = yaml.safe_load(ymlfile)

    # initialize droplets DataFrame from images
    droplets, rotation_theta = kchip_analyze.initialize_droplets(config)
    # print ('Rotation (degrees): ', rotation_theta*180/np.pi)

    droplets.to_csv(snakemake.output[0])
    # droplets.to_csv(snakemake.output[1])
    ## identify droplets in the same well from fit to masks
    # droplets = kchip_analyze.fit_droplets_to_mask(config,droplets,rotation_theta)
    # droplets.to_csv('droplets_assigned.csv')

    # remove overlap
    droplets = remove_overlap(config,droplets,show=0)
    unique_images = {tuple(row) for row in droplets[['IndexX','IndexY']].values}

    droplets['Well_ID']=0
    for x,y in unique_images:
        d = droplets.loc[(droplets['IndexX']==x) & (droplets['IndexY']==y)]
        positions = d[['ImageX','ImageY']].values
        same_well = droplets_in_same_well(positions)

        well_id = np.zeros((d.shape[0],1))

        for i,dab in enumerate(same_well):
            d_a, d_b = dab
            well_id[d_a]=i
            well_id[d_b]=i

        droplets.loc[(droplets['IndexX']==x) & (droplets['IndexY']==y),'Well_ID']=well_id

    # remove well_id = 0
    droplets = droplets[droplets['Well_ID']!=0]

    # apply hash
    droplets['Hash'] = droplets.apply(lambda row: hash((row['IndexX'],row['IndexY'],row['Well_ID'])),axis=1)

    # saves droplets assigned to wells and prints droplet count
    droplets.to_csv(snakemake.output[1])

    # removes wells with only one droplet and saves wells with combos
    droplets = droplets[droplets.groupby('Hash').Hash.transform(len) > 1]

    droplets.to_csv(snakemake.output[2])

    # corrects dye bleed through
    droplets['Dye 1 corr']=droplets['Dye 1'].values-(0.02*droplets['Dye 2'].values)
    droplets['Dye 2 corr']=droplets['Dye 2'].values-(0.14*droplets['Dye 1'].values)
    # reset column names for clustering code
    droplets[['R','G','B']]=\
        droplets[['Dye 0', 'Dye 1 corr', 'Dye 2 corr']]

    droplets.to_csv(snakemake.output[3])
