'''
Script: abx_only.py
Last update: 2021 March, mzhu

Description:
    1. abx-media growth norm se
    2. plotting into pdf with droplet level data
    3. plotting onto single axis
Inputs:
    1. data/output/{DAY}/{chip}_normalized.csv
    2. Keys-*.xlsx
    3. qc_params.csv
Outputs (csv):
    1. data/output/{DAY}/{chip}_abxonly_droplet_data.csv
    2. data/output/{DAY}/{chip}_abxonly_summary_data.csv
Outputs (figures):
    1. figures/{DAY}/{chip}_abxonly_MIC_curves.png
    2. figures/{DAY}/{chip}_abxonly_MIC_curves.pdf
'''
##############################################################################
################################## IMPORTS ###################################
import yaml
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
import statistics as st
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import multiprocessing as mp
from pathos.multiprocessing import ProcessingPool as Pool
import os, time
import warnings
warnings.filterwarnings('ignore')
from matplotlib.backends.backend_pdf import PdfPages
from plot_param import parameters
parameters()

##############################################################################
################################# FUNCTIONS ##################################
def bootstrap_median(ar):
    ''' resample from all dataframe rows '''
    return np.median(np.random.choice(ar, len(ar)))

def boot(df,boot_val,parallel=False):
    ''' boot_val: column name of values to bootstrap sample '''
    med_agg = bootstrap_median if parallel else 'median'
    grouped = df[['Label_left',boot_val]].groupby(['Label_left'])
    return grouped.aggregate(med_agg)

def boot_parallel(arglist):
	''' wrapper function for multiprocessing '''
	return boot(*arglist)

def MIC_boot_(df):
    '''
	bootstrap function to get standard error for MIC data
	df: normalized df with abx names mapped
    '''
    abxs = df.abx.unique()

    MED = df[['Label_left','norm_growth','abx',\
	'abx_name','abx_conc','abx_name_conc']].groupby(['Label_left',\
	'abx','abx_name','abx_conc','abx_name_conc']).median()

    cores = mp.cpu_count()//int(POOL_NUM)
    pool = Pool(cores)
    boot_ = pool.map(boot_parallel,([df,'norm_growth',j] for j in range(1,1001)))
    pool.close()
    pool.join()
    pool.clear()

    booted = pd.concat(boot_, axis=1)
    #print(booted.index)
    #print(MED.index)
    MED['norm_growth_SE'] = np.std(booted.values,axis=1)
    MED = MED.reset_index()
    SE_dict = dict(zip(MED.Label_left, MED['norm_growth_SE']))
    df['norm_growth_SE'] = df['Label_left'].map(SE_dict)
    return df,MED

def plot_MIC(df,out_path,chip,t=0):
	''' plotting MIC data, each abx as a panel in multipage PDF '''
	abxs = df.abx_name.unique()
	pp = PdfPages(out_path+'_abxonly_MIC_curves.pdf')
	for i in abxs:
		abx_drops = df[df.abx_name==i]
		abx_drops.sort_values(by='abx_conc',inplace=True)
		med = abx_drops.groupby('abx_conc').median()
		count = abx_drops.groupby('abx_conc').count()
		plt.figure()
		sns.scatterplot(data=abx_drops, x="abx_conc", y='norm_growth', hue="abx_name_conc",alpha=0.3)
		plt.errorbar(med.index,med['norm_growth'],yerr=med['norm_growth_SE'],c='black',capsize=3)
		for j in np.sort(med.index):
        # for j in ['a','b','c']:
			plt.text(j,med.loc[j,'norm_growth']+0.4,str(count.loc[j,'norm_growth_SE']))
		plt.title(chip+' - '+i)
		plt.ylim([0,2.5])
		plt.legend(bbox_to_anchor=(1.2,1))
		plt.savefig(pp, format='pdf',bbox_inches='tight',dpi=300)
	pp.close()
################################################################################
qc_params = snakemake.input[2]
# get qc_params
qc_paramsdf = pd.read_csv(qc_params)
POOL_NUM = qc_paramsdf.iloc[4,1]

norm_file = snakemake.input[0]
keys = pd.read_excel(snakemake.input[1],index_col=0,sheet_name='ABX')

norm = pd.read_csv(norm_file,index_col=0)
day = norm_file.split('/')[-1].split('_')[0]
chip = norm_file.split('/')[-1].split('_normalized')[0]
data_folder = 'data/output/'+day+'/'+chip
figure_folder='figures/'+day+'/'+chip

# create strain-abx dict from keys
ab1_key = dict(zip(keys.index, keys.ab1))
ab2_key = dict(zip(keys.index, keys.ab2))
# ab3_key = dict(zip(keys.index, keys.ab3))
# ab4_key = dict(zip(keys.index, keys.ab4))
kp1_key = dict(zip(keys.index, keys.kp1))
kp2_key = dict(zip(keys.index, keys.kp2))
# kp3_key = dict(zip(keys.index, keys.kp3))
# kp4_key = dict(zip(keys.index, keys.kp4))
pa1_key = dict(zip(keys.index, keys.pa1))
pa2_key = dict(zip(keys.index, keys.pa2))
# pa3_key = dict(zip(keys.index, keys.pa3))
# pa4_key = dict(zip(keys.index, keys.pa4))
strains = ['ab1','ab2','kp1','kp2','pa1','pa2']
# strains = ['ab1','ab2','ab3','ab4','kp1','kp2','kp3','kp4','pa1','pa2','pa3','pa4']
strain_key = dict(zip(strains,[ab1_key,ab2_key,kp1_key,kp2_key,pa1_key,pa2_key]))
# strain_key = dict(zip(strains,[ab1_key,ab2_key,ab3_key,ab4_key,kp1_key,kp2_key,kp3_key,kp4_key,pa1_key,pa2_key,pa3_key,pa4_key]))

STRAIN = norm_file.split('/')[-1].split('_')[1]

# map abx names
abxdict = strain_key[STRAIN]
df = norm[(norm.Label_left.str.contains('abx'))&(norm.Label_right.str.contains('CAMHB'))]
df['abx'] = df.Label_left.apply(lambda x: x[:-1])
df['abx_name'] = df['abx'].map(abxdict)
df['abx_conc'] = df.Label_left.apply(lambda x: x[-1])
df['abx_name_conc'] = df.abx_name+'_'+df.abx_conc

# bootstrap for se
abxnorm,abxMED = MIC_boot_(df)

# sort based on abx names and conc
abxnorm.sort_values(by='abx_name_conc',inplace=True)
abxMED.sort_values(by='abx_name_conc',inplace=True)

# add chip info
abxnorm['chipID'] = chip
abxMED['chipID'] = chip

# save output
abxnorm.to_csv(snakemake.output[0])
abxMED.to_csv(snakemake.output[1])

# plot abx into pdf
plot_MIC(abxnorm,figure_folder,chip)

# plot all abx on single axeis
num = abxMED.abx_name.unique().shape[0]
clr_map = plt.cm.viridis # LinearSegmentedColormap
N_clrs = min(clr_map.N,num)
map_clrs = [clr_map(int(x*clr_map.N/N_clrs)) for x in range(N_clrs)]
plt.figure(figsize=(20,3))
# for i in np.arange(0, num, 3):
#     plt.errorbar(x=abxMED.abx_name_conc[i:i+3], y=abxMED.norm_growth[i:i+3], \
# 	yerr=abxMED.norm_growth_SE[i:i+3], c=map_clrs[i])
# changed
xticks_ = [0]
val = 0
for i in enumerate(abxMED.abx_name.unique()):
    plt.errorbar(x=abxMED.loc[(abxMED.abx_name==i[1]),'abx_name_conc'], y=abxMED.loc[(abxMED.abx_name==i[1]),'norm_growth'], \
	yerr=abxMED.loc[(abxMED.abx_name==i[1]),'norm_growth_SE'], c=map_clrs[i[0]])
    val+=abxMED[abxMED.abx_name==i[1]].shape[0]
    xticks_.append(val)

unique = np.unique(abxMED.abx_name)
# plt.xticks(np.arange(0, num, 3), unique, rotation=45)
plt.xticks(xticks_[:-1], unique, rotation=45)
plt.ylabel('median relative growth')
plt.xlabel('antibiotics in descending concentration order')
plt.title(chip+' - all abx curves')
plt.savefig(figure_folder+'_abxonly_MIC_curves.png')
