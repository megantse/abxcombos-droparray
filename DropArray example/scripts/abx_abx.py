'''
Script: abx_abx.py
Last update: 2021 March, mzhu

Description:
    1. abx-abx bliss score and se
    2. abx-abx null bliss score and se
    3. t-test for signficance
    4. synergy control bliss score and se (pval not calc)
Inputs:
    1. data/output/{DAY}/{chip}_normalized.csv
    2. Keys-*.xlsx
    3. synergy_control.xlsx
    4. qc_params.csv
Outputs (csv):
    1. data/output/{DAY}/{chip}_abxabx_droplet_data.csv
    2. data/output/{DAY}/{chip}_abxabx_summary_data.csv
    3. data/output/{DAY}/{chip}_abxabx_negctrl_droplet_data.csv
    4. data/output/{DAY}/{chip}_abxabx_negctrl_summary_data.csv
    5. data/output/{DAY}/{chip}_abxabx_negctrl_summed_summary_data.csv
    6. data/output/{DAY}/{chip}_abxabx_summed_summary_data.csv
    7. data/output/{DAY}/{chip}_abxabx_synctrl.csv
Outputs (figures):
    1. figures/{DAY}/{chip}_abxabx_bliss_distribution.png # sums
    2. figures/{DAY}/{chip}_abxabx_bliss_volcano.png # sums
    3. figures/{DAY}/{chip}_abxabx_blissSE_heatmap.png # sums
    4. figures/{DAY}/{chip}_abxabx_normgrowth_heatmap.png # individual
    5. figures/{DAY}/{chip}_abxabx_null_distribution.png # sums (x3)
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
def separate_sorted_df(df):
    '''
    Sorts normalized df into combo types
    df: normalized, abx cp name mapping not needed
    '''
    sep_dict = dict()
    abx = df[df['Label_left'].str.contains('abx')]
    sep_dict['abx_pair'] = abx[abx['Label_right'].str.contains('abx')] # abx-abx
    sep_dict['abx_bugs'] = abx[abx['Label_right'].str.contains('BUGS')]
    sep_dict['abx_media'] = abx[abx['Label_right'].str.contains('CAMHB')]
    sep_dict['abx_cp'] = abx[abx['Label_right'].str.contains('cp')]

    bugs = df[df['Label_left'].str.contains('BUGS')]
    sep_dict['bug_pair'] = bugs[bugs['Label_right'].str.contains('BUGS')] # bug-bug
    sep_dict['bug_media'] = bugs[bugs['Label_right'].str.contains('CAMHB')]
    sep_dict['bug_cp'] = bugs[bugs['Label_right'].str.contains('cp')]

    media = df[df['Label_left'].str.contains('CAMHB')]
    sep_dict['media_pair'] = media[media['Label_right'].str.contains('CAMHB')]

    return sep_dict

def bootstrap_median(df):
    ''' Resample from all dataframe rows '''
    return np.median(np.random.choice(df, len(df)))

def bliss_calc(abx_data, ctrl_a, ctrl_b, parallel=False, save=False):
    '''
    Calculates Bliss scores for given abx abx pair.

    Inputs:
    :abx_data: (DataFrame) abx-abx combo data, normalized to bug-bug
    :ctrl_a: (DataFrame) abx-bug combo data, normalized to bug-bug
    :ctrl_b: (DataFrame) abx-bug combo data, normalized to bug-bug

    Outputs:
    if parallel: just the bliss column returned
    if non-parallel: entire bliss df returned
    '''

    med_agg = bootstrap_median if parallel else 'median'

    abg = ctrl_a[['Label_left','norm_growth_inh']].groupby('Label_left')
    control_a = abg.aggregate(med_agg) # controls grouped by abx

    ctrl_b = ctrl_b.rename(columns={'Label_right':'BUGS', 'Label_left':'Label_right'})
    cbg = ctrl_b[['Label_right','norm_growth_inh']].groupby('Label_right')
    control_b = cbg.aggregate(med_agg) # controls grouped by abx

    data_df = abx_data[['Label_left','Label_right','norm_growth_inh']].groupby(['Label_left','Label_right'])
    data_df = data_df.aggregate(med_agg)

    synergy_df = data_df.merge(control_a,left_index=True,right_index=True)
    synergy_df.rename(columns={'norm_growth_inh_y':'E_a'}, inplace=True)

    synergy_df = synergy_df.merge(control_b,left_index=True,right_index=True)
    synergy_df.rename(columns={'norm_growth_inh':'E_b', 'norm_growth_inh_x': 'Eab'}, inplace=True)
    bliss = synergy_df.copy()
    bliss['bliss'] = (bliss.Eab - (bliss.E_a + bliss.E_b - bliss.E_a*bliss.E_b)).astype('float')

    if parallel:
        return bliss.loc[:, 'bliss']
    else:
        return bliss

def bliss_parallel(arglist):
    return bliss_calc(*arglist)

##############################################################################
# read in normalized, keys, and synergy ctrl files and set other vars
norm_file = snakemake.input[0]
keys = pd.read_excel(snakemake.input[1],index_col=0,sheet_name='ABX')
syn_ctrl = pd.read_excel(snakemake.input[2],index_col=0)

qc_params = snakemake.input[3]
# get qc_params
qc_paramsdf = pd.read_csv(qc_params)
POOL_NUM = qc_paramsdf.iloc[4,1]

norm = pd.read_csv(norm_file,index_col=0)
day = norm_file.split('/')[-1].split('_')[0]
chip = norm_file.split('/')[-1].split('_normalized')[0]
data_folder = 'data/output/'+day+'/'+chip
figure_folder='figures/'+day+'/'+chip

STRAIN = norm_file.split('/')[-1].split('_')[1]

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


##############################################################################
############################## ABX ABX BLISS #################################
# create dictionary of df combo types
bliss_dict = separate_sorted_df(norm)

# parallelized bliss scoring
# create the multiprocessing pool
cores = mp.cpu_count()//int(POOL_NUM)
pool = Pool(cores)

# process the dataFrame by mapping function to each df across the pool
bliss_pool = pool.map(bliss_parallel,
                    ([bliss_dict['abx_pair'],
                    bliss_dict['abx_bugs'],
                    bliss_dict['abx_bugs'],
                    i] for i in range(1,1001)))

# close down the pool and join
pool.close()
pool.join()
pool.clear()

# calc bliss and se
bliss_boot = pd.concat(bliss_pool, axis=1)
med_bliss = bliss_boot.transpose().median().reset_index()
med_bliss['comboID'] = list(zip(med_bliss.Label_left,med_bliss.Label_right))
med_dict = dict(zip(med_bliss.comboID,med_bliss.iloc[:,2]))

se_bliss = bliss_boot.transpose().std().reset_index()
se_bliss['comboID'] = list(zip(se_bliss.Label_left,se_bliss.Label_right))
se_dict = dict(zip(se_bliss.comboID,se_bliss.iloc[:,2]))

# add bliss and se to abx abx df
abxabx = bliss_dict['abx_pair']
abxabx['comboID'] = list(zip(abxabx.Label_left,abxabx.Label_right))
abxabx['bliss_med'] = abxabx['comboID'].map(med_dict)
abxabx['bliss_se'] = abxabx['comboID'].map(se_dict)

# map abx names
abxdict = strain_key[STRAIN]
abxabx['abx1'] = abxabx.Label_left.apply(lambda x: x[:-1])
abxabx['abx2'] = abxabx.Label_right.apply(lambda x: x[:-1])
abxabx['abx_name1'] = abxabx['abx1'].map(abxdict)
abxabx['abx_name2'] = abxabx['abx2'].map(abxdict)
abxabx['abx_conc1'] = abxabx.Label_left.apply(lambda x: x[-1])
abxabx['abx_conc2'] = abxabx.Label_right.apply(lambda x: x[-1])
abxabx['abx_name_conc1'] = abxabx.abx_name1+'_'+abxabx.abx_conc1
abxabx['abx_name_conc2'] = abxabx.abx_name2+'_'+abxabx.abx_conc2

# calc tstat
abxabx['bliss_tstat'] = abxabx.bliss_med / abxabx.bliss_se

# save droplet data
abxabx['chipID'] = chip
abxabx.to_csv(snakemake.output[0])

# calc and save summary data
abxabx_med = abxabx[['abx_name1','abx_name2','abx_name_conc1','abx_name_conc2',\
'bliss_med','bliss_se','bliss_tstat']].groupby(['abx_name1','abx_name2',\
'abx_name_conc1','abx_name_conc2']).median().reset_index()
abxabx_med['chipID'] = chip
abxabx_med.to_csv(snakemake.output[1])

# calc summed bliss scores
abxabx_sum = abxabx_med[['abx_name1','abx_name2','bliss_med',\
'bliss_se']].groupby(['abx_name1','abx_name2']).sum().reset_index()
abxabx_sum['bliss_tstat'] = abxabx_sum.bliss_med / abxabx_sum.bliss_se

# plotting summed bliss score distribution
plt.figure()
plt.title(chip+' - abx-abx bliss sum distribution')
sns.distplot(abxabx_sum['bliss_med'])
plt.savefig(figure_folder+'_abxabx_bliss_distribution.png')

##############################################################################
######################## PLOTTING HEATMAPS ###################################
def matrixize(df,ind,col,val):
    '''
    Makes melted df into a filled square matrix
    Inputs:
    :df: dataframe in melted form
    :ind: column name to make vals --> index
    :col: column name to make vals --> columns
    :val: values to fill index x columns 2D matrix
    Output:
    :df_mat: matrixcized df for heatmap plotting
    '''
    df_pivot = pd.DataFrame(index=df[ind].unique(),columns=df[col].unique())
    for i in df_pivot.index:
        for j in df_pivot.columns:
            try:
                df_pivot.loc[i,j] = df[(df[ind]==i)&(df[col]==j)][val].values[0]
                df_pivot.loc[j,i] = df[(df[ind]==i)&(df[col]==j)][val].values[0]
            except:
                pass
    df_pivot = df_pivot.fillna(0)
    df_pivot.sort_index(inplace=True)
    df_pivot = df_pivot.reindex(sorted(df_pivot.columns), axis=1)
    return df_pivot

def matrix_mask(df):
    '''
    Creates mask for heatmap of square matrix
    :df: matrixized dataframe
    :mask: mask to avoid double plotting values in a square heatmap
    '''
    mask = np.zeros_like(df, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True
    return mask

# plotting summed bliss score heatmaps
df_bliss = matrixize(abxabx_sum,ind='abx_name1',col='abx_name2',val='bliss_med')
mask_bliss = matrix_mask(df_bliss)
df_se = matrixize(abxabx_sum,ind='abx_name1',col='abx_name2',val='bliss_se')
mask_se = matrix_mask(df_se)

fig, ax = plt.subplots(1,2,figsize=(20,8))
ax[0].set_title(chip+' - abx-abx summed bliss scores')
ax[1].set_title(chip+' - abx-abx summed bliss SE')
sns.heatmap(df_bliss,cmap='coolwarm',mask=mask_bliss,cbar_kws={'pad':0.01,'aspect':50},ax=ax[0])
sns.heatmap(df_se,cmap='coolwarm',mask=mask_se,cbar_kws={'pad':0.01,'aspect':50},ax=ax[1])
plt.savefig(figure_folder+'_abxabx_blissSE_heatmap.png')

# plotting growth heatmaps
abxabx_growth = abxabx[['abx_name_conc1','abx_name_conc2',\
'norm_growth']].groupby(['abx_name_conc1','abx_name_conc2']).median().reset_index()
abxabx_ct = abxabx[['abx_name_conc1','abx_name_conc2',\
'norm_growth']].groupby(['abx_name_conc1','abx_name_conc2']).count().reset_index()

abxabx_growth_ = matrixize(abxabx_growth,'abx_name_conc1','abx_name_conc2','norm_growth')
mask_inh = matrix_mask(abxabx_growth_)

abxabx_ct_ = matrixize(abxabx_ct,'abx_name_conc1','abx_name_conc2','norm_growth')
# abxabx_ct_.to_csv('abxabx_count.csv')
mask_ct = matrix_mask(abxabx_ct_)
# print(mask_ct)

fig, ax = plt.subplots(1,2,figsize=(40,16))
ax[0].set_title(chip+' - abx-abx normalized growth')
ax[1].set_title(chip+' - abx-abx counts')
sns.heatmap(abxabx_growth_,cmap='coolwarm',mask=mask_inh,cbar_kws={'pad':0.01,'aspect':50},ax=ax[0])
sns.heatmap(abxabx_ct_.astype('float'),cmap='coolwarm',mask=mask_ct,annot=True,cbar_kws={'pad':0.01,'aspect':50},ax=ax[1])
plt.savefig(figure_folder+'_abxabx_normgrowth_heatmap.png')

##############################################################################
######################## NEGATIVE CONTROLS ###################################
def ctrl_bliss_calc(abx_bug1, abx_bug2, bug_pair, parallel=False, save=False):
    '''
    Calculates Bliss scores for abx abx null distribution

    Inputs:
    :abx_bug1: (DataFrame) abx-bug combo data, normalized to bug-bug
    :abx_bug2: (DataFrame) abx-bug combo data, normalized to bug-bug
    :bug_pair: (DataFrame) bug-bug combo data, normalized to bug-bug

    Outputs:
    if parallel: just the bliss scores
    if non-parallel: entire DataFrame with E_a, E_b, and bliss columns added
    '''

    med_agg = bootstrap_median if parallel else 'median'

    abg = abx_bug1[['Label_left','norm_growth_inh']].groupby('Label_left')
    control_a = abg.aggregate(med_agg) # controls grouped by abx

    cbg = bug_pair.norm_growth_inh.values
    control_b = np.median(np.random.choice(cbg,len(cbg)))

    data_df = abx_bug2[['Label_left','norm_growth_inh']].groupby(['Label_left'])
    data_df = data_df.aggregate(med_agg)

    synergy_df = data_df.merge(control_a,left_index=True,right_index=True)
    synergy_df.rename(columns={'norm_growth_inh_x': 'Eab','norm_growth_inh_y':'E_a'}, inplace=True)

    synergy_df['E_b'] = control_b
    bliss = synergy_df.copy()
    bliss['bliss'] = (bliss.Eab - (bliss.E_a + bliss.E_b - bliss.E_a*bliss.E_b)).astype('float')

    if parallel:
        return bliss.loc[:, 'bliss']
    else:
        return bliss

def ctrl_bliss_parallel(arglist):
    return ctrl_bliss_calc(*arglist)

# extract negative control pairs
neg_Eab = bliss_dict['abx_bugs'].sort_values(by=['Label_left','Label_right']).iloc[::2,:]
neg_Ea = bliss_dict['abx_bugs'].sort_values(by=['Label_left','Label_right']).iloc[1::2,:]
neg_Eb = bliss_dict['bug_pair']

# parallelized Bliss scoring
# create the multiprocessing pool
cores = mp.cpu_count()//int(POOL_NUM)
pool = Pool(cores)

# process the DataFrame by mapping function to each df across the pool
ctrl_bliss_pool = pool.map(ctrl_bliss_parallel,([neg_Eab,neg_Ea,neg_Eb,i] for i in range(1,1001)))

# close down the pool and join
pool.close()
pool.join()
pool.clear()

# calc bliss and se for negative controls
ctrl_bliss_boot = pd.concat(ctrl_bliss_pool, axis=1)
ctrl_med_bliss = ctrl_bliss_boot.transpose().median().reset_index()
ctrl_med_dict = dict(zip(ctrl_med_bliss.Label_left,ctrl_med_bliss.iloc[:,1]))

ctrl_se_bliss = ctrl_bliss_boot.transpose().std().reset_index()
ctrl_se_dict = dict(zip(ctrl_se_bliss.Label_left,ctrl_se_bliss.iloc[:,1]))

# add bliss and se to abx - bugs df
abxbug = bliss_dict['abx_bugs']
abxbug['bliss_med'] = abxbug['Label_left'].map(ctrl_med_dict)
abxbug['bliss_se'] = abxbug['Label_left'].map(ctrl_se_dict)

# map abx names
abxdict = strain_key[STRAIN]
abxbug['abx'] = abxbug.Label_left.apply(lambda x: x[:-1])
abxbug['abx_name'] = abxbug['abx'].map(abxdict)
abxbug['abx_conc'] = abxbug.Label_left.apply(lambda x: x[-1])
abxbug['abx_name_conc'] = abxbug.abx_name+'_'+abxbug.abx_conc

# calc tstat
abxbug['bliss_tstat'] = abxbug.bliss_med / abxbug.bliss_se

# calc median bliss score and SE
abxbug_med = abxbug[['abx_name_conc','bliss_med',\
'bliss_se','bliss_tstat']].groupby('abx_name_conc').median().reset_index()
abxbug_med['abx_name'] = abxbug_med.abx_name_conc.apply(lambda x: x.split('_')[0])

# calc 3 sum bliss scores and SE
abxbug_sum = abxbug_med[['abx_name','bliss_med',\
'bliss_se']].groupby('abx_name').sum()

# calc summed tstat
abxbug_sum['bliss_tstat'] = abxbug_sum.bliss_med / abxbug_sum.bliss_se

# multiple null distribution by factor of 3
abxbug_sum['bliss_tstat_x3'] = abxbug_sum.bliss_tstat*3

# add chip info
abxbug['chipID'] = chip
abxbug_med['chipID'] = chip
abxbug_sum['chipID'] = chip

# save outputs
abxbug.to_csv(snakemake.output[2])
abxbug_med.to_csv(snakemake.output[3])
abxbug_sum.to_csv(snakemake.output[4])

# fit null to t distribution
null = abxbug_sum['bliss_tstat_x3'].values.tolist()
tdof, tloc, tscale = stats.t.fit(null, floc = 0.)
null_fit = stats.t.rvs(tdof, loc=tloc, scale=tscale, size=10000)
rv = stats.t(tdof, loc=tloc, scale=tscale)

# plotting tstat distribution
a = np.arange(-5,5,0.1)
plt.figure()
plt.title(chip+' - abx-abx null distribution')
plt.plot(a, rv.pdf(a), '-', label='fitted data')
plt.hist(abxbug_sum['bliss_tstat_x3'],density=True,label='empirical data')
plt.legend(bbox_to_anchor=(1.25,1))
plt.savefig(figure_folder+'_abxabx_null_distribution.png')

def t_test(tstat):
    ''' Get p value from t test '''
    global tdof, tloc, tscale
    return stats.t.sf(abs(tstat), tdof, loc=tloc, scale=tscale)

# calc p values from t-test
abxabx_sum['bliss_pval'] = abxabx_sum['bliss_tstat'].apply(t_test)
abxabx_sum['bliss_log10pval'] = -np.log10(abxabx_sum['bliss_pval'])

# save summed bliss data
abxabx_sum['chipID'] = chip
abxabx_sum.to_csv(snakemake.output[5])

# plot volcano for abx - abx pairs
plt.figure()
sns.scatterplot(data=abxabx_sum,x='bliss_med',y='bliss_log10pval',alpha=0.7)
plt.title(chip+' - abx-abx volcano plot')
plt.savefig(figure_folder+'_abxabx_bliss_volcano.png')

##############################################################################
######################## SYNERGY CONTROLS ###################################
# make synergy dictionary
synctrl_dict = syn_ctrl.to_dict('index')
strain_ctrl = synctrl_dict[STRAIN]

# make sum type dictionary
sumtype = {syn_ctrl.loc[STRAIN,'abx1_name']:syn_ctrl.loc[STRAIN,'abx1_conc'],\
syn_ctrl.loc[STRAIN,'abx2_name']:syn_ctrl.loc[STRAIN,'abx2_conc']}

# create abx and conc columns
abxabx_med['abx1_name'] = abxabx_med['abx_name_conc1'].apply(lambda x: x.split('_')[0])
abxabx_med['abx2_name'] = abxabx_med['abx_name_conc2'].apply(lambda x: x.split('_')[0])
abxabx_med['abx1_conc'] = abxabx_med['abx_name_conc1'].apply(lambda x: x.split('_')[-1])
abxabx_med['abx2_conc'] = abxabx_med['abx_name_conc2'].apply(lambda x: x.split('_')[-1])

# pick out control pairs
df_R = abxabx_med[(abxabx_med.abx1_name==strain_ctrl['abx1_name'])&(abxabx_med.abx2_name==\
strain_ctrl['abx2_name'])]
df_L = abxabx_med[(abxabx_med.abx2_name==strain_ctrl['abx1_name'])&(abxabx_med.abx1_name==\
strain_ctrl['abx2_name'])]
synctrl_df = pd.concat([df_R,df_L])

# redefine abx1 and abx2
abx1 = synctrl_df.abx1_name.values[0]
abx2 = synctrl_df.abx2_name.values[0]

# get sum type
if sumtype[abx1]!= 'all':
    synctrl_sum = synctrl_df.loc[(synctrl_df.abx1_conc==sumtype[abx1]),['bliss_med',\
    'bliss_se']].sum().reset_index().T
    synctrl_sum.columns = synctrl_sum.iloc[0]
    synctrl_sum.drop(synctrl_sum.index[0], inplace = True)
elif sumtype[abx2]!= 'all':
    synctrl_sum = synctrl_df.loc[(synctrl_df.abx2_conc==sumtype[abx2]),['bliss_med',\
    'bliss_se']].sum().reset_index().T
    synctrl_sum.columns = synctrl_sum.iloc[0]
    synctrl_sum.drop(synctrl_sum.index[0], inplace = True)

# create abx and conc columns
synctrl_sum['abx1_name'] = abx1
synctrl_sum['abx1_conc'] = sumtype[abx1]
synctrl_sum['abx2_name'] = abx2
synctrl_sum['abx2_conc'] = sumtype[abx2]
synctrl_sum['bliss_tstat'] = synctrl_sum.bliss_med / synctrl_sum.bliss_se
synctrl_sum['chipID'] = chip
synctrl_sum.to_csv(snakemake.output[6])
