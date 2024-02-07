'''
Script: cp_only.py
Last update: 2024 Feb, mtse

Description:
    1. bug-cp growth norm se
    2. pick out high activity cp
    3. plot growth (inh) as histogram
Inputs:
    1. data/output/{DAY}/{chip}_normalized.csv
    2. Keys-*.xlsx
    3. qc_params.csv
Outputs (csv):
    1. data/output/{DAY}/{chip}_cponly_droplet_data.csv
    2. data/output/{DAY}/{chip}_cponly_summary_data.csv
    3. data/output/{DAY}/{chip}_cponly_high_activity.csv
Outputs (figures):
    1. figures/{DAY}/{chip}_cponly_histogram_normgrowth.png
    2. figures/{DAY}/{chip}_cponly_histogram_normgrowthinh.png
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
    ''' Resample from all dataframe rows '''
    return np.median(np.random.choice(ar, len(ar)))

def boot(df,boot_val,parallel=False):
    '''
    :df: normalized df
    :boot_val: column name of values to bootstrap sample '''
    med_agg = bootstrap_median if parallel else 'median'
    grouped = df[['Label_right',boot_val]].groupby(['Label_right'])
    return grouped.aggregate(med_agg)

def boot_parallel(arglist):
	''' wrapper function for multiprocessing '''
	return boot(*arglist)

def CP_boot_(df):
    '''
	Bootstrap function to get standard error for cp only data
	df: normalized df with cp names mapped
    '''

    MED = df[['Label_right','norm_growth','cp_name']].groupby(['Label_right','cp_name']).median()
    cores = mp.cpu_count()//int(POOL_NUM)
    pool = Pool(cores)
    boot_ = pool.map(boot_parallel,([df,'norm_growth',j] for j in range(1,1001)))
    pool.close()
    pool.join()
    pool.clear()

    booted = pd.concat(boot_, axis=1)
    MED['norm_growth_SE'] = np.std(booted.values,axis=1)
    MED = MED.reset_index()
    SE_dict = dict(zip(MED.Label_right, MED['norm_growth_SE']))
    df['norm_growth_SE'] = df['Label_right'].map(SE_dict)
    return df,MED

def t_pval(df, med_bugs_df):
    '''
    Calculate the p-value of samples growth norm values being less than the
    control, bug-media
    Function used for groupby apply
    df: normalized df with 1X bugs
    med_bugs_df: normalized df with bug-media pairs
    '''

    sample = df['norm_growth']
    null = med_bugs_df['norm_growth']
    _,p_val = stats.ttest_ind(sample,null,alternative='less')
    return p_val
################################################################################
norm_file = snakemake.input[0]
keys = pd.read_excel(snakemake.input[1],index_col=0,sheet_name='CP')

qc_params = snakemake.input[2]
# get qc_params
qc_paramsdf = pd.read_csv(qc_params)
POOL_NUM = qc_paramsdf.iloc[4,1]

norm = pd.read_csv(norm_file,index_col=0)
day = norm_file.split('/')[-1].split('_')[0]
chip = norm_file.split('/')[-1].split('_normalized')[0]
data_folder = 'data/output/'+day+'/'+chip
figure_folder='figures/'+day+'/'+chip

STRAIN = norm_file.split('/')[-1].split('_')[1]
CPP = norm_file.split('/')[-1].split('_')[2]

# map cp names
df_ = norm[(norm.Label_left.str.contains('BUGS'))&(norm.Label_right.str.contains('cp'))]
controls = norm[(norm.Label_left.str.contains('CAMHB'))&(norm.Label_right.str.contains('CAMHB'))]
df = pd.concat([controls,df_])
cppkeys = keys[keys.Blainey_cpp==CPP]
cppdict = dict(zip(cppkeys.Blainey_cp,cppkeys.Broad_Sample))
df['cp_name'] = df['Label_right'].map(cppdict)
df.loc[(df.Label_left.str.contains('CAMHB'))&(df.Label_right.str.contains('CAMHB')),'cp_name'] = 'CAMHB'

# bootstrap for norm growth se
cpnorm,cpMED = CP_boot_(df)
cpMED['norm_growth_inh'] = 1 - cpMED.norm_growth

# sort based on cp name
cpnorm.sort_values(by='Label_right',inplace=True)
cpMED.sort_values(by='Label_right',inplace=True)

# add chip info
cpnorm['chipID'] = chip
cpMED['chipID'] = chip

# plot growth (inh) histograms
f, (ax_hist, ax_box) = plt.subplots(2, sharex=True,
                                    gridspec_kw={"height_ratios": (.85, .15)})

sns.distplot(cpMED.norm_growth, ax=ax_hist)
sns.boxplot(cpMED.norm_growth, ax=ax_box)
ax_box.set(yticks=[])
sns.despine(ax=ax_hist)
sns.despine(ax=ax_box, left=True)
plt.suptitle(chip+' - compound only activity')
plt.savefig(figure_folder+'_cponly_histogram_normgrowth.png')

f, (ax_hist, ax_box) = plt.subplots(2, sharex=True,
                                    gridspec_kw={"height_ratios": (.85, .15)})

sns.distplot(cpMED.norm_growth_inh, ax=ax_hist)
sns.boxplot(cpMED.norm_growth_inh, ax=ax_box)
ax_box.set(yticks=[])
sns.despine(ax=ax_hist)
sns.despine(ax=ax_box, left=True)
plt.suptitle(chip+' - compound only activity')
plt.savefig(figure_folder+'_cponly_histogram_normgrowthinh.png')

# calculate p-values and place into cpMED dfs
bugmed_df = norm[(norm.Label_left.str.contains('BUGS'))&(norm.Label_right.str.contains('CAMHB'))]
cp_pval = cpnorm.groupby('Label_right').apply(t_pval,med_bugs_df=bugmed_df).reset_index()
cpMED = cpMED.merge(cp_pval,on=['Label_right']).rename(columns={0:'normgrowth_pval'})
cpMED['normgrowth_log10pval'] = -np.log10(cpMED['normgrowth_pval'])

plt.figure()
sns.scatterplot(data=cpMED, x='norm_growth',y='normgrowth_log10pval')
plt.title(chip + ' - compound only activity compared to controls')
plt.xlabel('relative growth')
plt.ylabel('-log10 pval')
plt.savefig(figure_folder+'_cponly_scatterplot_normgrowth.png')

# save output
cpnorm.to_csv(snakemake.output[0])
cpMED.to_csv(snakemake.output[1])

# pick out cp with growth inh > 30%
CPhigh_act = cpMED.loc[cpMED.norm_growth_inh>=30,:]
CPhigh_act.to_csv(snakemake.output[2])
