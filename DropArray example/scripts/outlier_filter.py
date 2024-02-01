'''
Script: outlier_filter.py
Last update: 2021 March, mzhu

Description: OUTLIER FILTERING for every combo type
Input:
    1. data/output/{DAY}/{chip}_qcfiltered.csv
Output:
    1. data/output/{DAY}/{chip}_outlierfiltered.csv
'''
##############################################################################
################################## IMPORTS ###################################
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import multiprocessing as mp
from pathos.multiprocessing import ProcessingPool as Pool
import warnings
warnings.filterwarnings('ignore')

##############################################################################
########################### OUTLIER FILTERING ###############################
qcfil_file = snakemake.input[0]

qcfil = pd.read_csv(qcfil_file,index_col=0)
day = qcfil_file.split('/')[-1].split('_')[0]
chip = qcfil_file.split('/')[-1].split('_qcfiltered')[0]
data_folder = 'data/output/'+day+'/'+chip
figure_folder='figures/'+day+'/'+chip

# x1,y1 = qcfil.shape;print(x1,y1)
IQR_df = qcfil.groupby(['Label_left','Label_right'])['t0_norm2'].quantile([0.25, 0.75]).reset_index()#.unstack(level=1)
IQR_df['comboID'] = list(zip(IQR_df.Label_left,IQR_df.Label_right))
Q1 = dict(zip(IQR_df.loc[(IQR_df.level_2==0.25),'comboID'],IQR_df.loc[(IQR_df.level_2==0.25),'t0_norm2']))
Q3 = dict(zip(IQR_df.loc[(IQR_df.level_2==0.75),'comboID'],IQR_df.loc[(IQR_df.level_2==0.75),'t0_norm2']))

qcfil['comboID'] = list(zip(qcfil.Label_left,qcfil.Label_right))
qcfil['Q1'] = qcfil['comboID'].map(Q1)
qcfil['Q3'] = qcfil['comboID'].map(Q3)

qcfil['IQR'] = qcfil.Q3 - qcfil.Q1
qcfil['lower_range'] = qcfil.Q1 - (1.5*qcfil.IQR)
qcfil['upper_range'] = qcfil.Q3 + (1.5*qcfil.IQR)
qcfil = qcfil[(qcfil.t0_norm2>=qcfil.lower_range)&(qcfil.t0_norm2<=qcfil.upper_range)]
# x2,y2 = qcfil.shape;print(x2,y2)
# print((x1 - x2)/x1*100,'percent filtered from outliers')
qcfil['chipID'] = chip

# save outputs
qcfil.to_csv(snakemake.output[0])
