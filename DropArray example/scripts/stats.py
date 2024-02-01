'''
Script: stats.py
Last update: 2021 March, mzhu

Description: STATS
    PDF 1 - # droplets and microwells, cluster counts,
            control distribution and counts, correlation, standard z primes,
            abx z primes, # microwells per combo type, droplet location,
            compound activity, abx curves, abx-abx bliss, abx-cp growth
    PDF 2 - chip hits, # droplets and microwells, cluster counts,
            control distribution and counts, correlation, standard z primes,
            abx z primes, # microwells per combo type, droplet location,
            compound activity, abx curves, abx-abx bliss, abx-cp growth
    PDF 3 - all PNG's
Input:
    1. data/output/{DAY}/{chip}_normalized.csv
    2. 'Keys-secondary.xlsx'
    3. data/interim/{DAY}/{chip}_missing_clusters.csv
    4. data/interim/{DAY}/{chip}_droplets_b.csv
    5. data/interim/{DAY}/{chip}_clustered.csv
    6. data/interim/{DAY}/{chip}_trimmed.csv
    7. data/output/{DAY}/{chip}_qcfiltered.csv
    8. data/output/{DAY}/{chip}_outlierfiltered.csv
    9. data/interim/{DAY}/{chip}_pre_post.csv
    10. data/output/{DAY}/{chip}_cponly_summary_data.csv
    11. data/output/{DAY}/{chip}_abxonly_summary_data.csv
    12. data/output/{DAY}/{chip}_abxabx_summed_summary_data.csv
Output (csv): # missing file is re-written as well
    1. data/output/{DAY}/{chip}_stat_zprime_st.csv
	2. data/output/{DAY}/{chip}_stat_zprime_abx.csv
Output (figures):
    1. figures/{DAY}/{chip}_stats_compiled.pdf
	2. figures/{DAY}/{chip}_compiled_chip_hitstats.pdf
    3. figures/{DAY}/{chip}_compiled_all_figures.pdf
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
import glob
import multiprocessing as mp
from pathos.multiprocessing import ProcessingPool as Pool
import os, time
from PIL import Image
from PyPDF2 import PdfFileMerger
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')
from matplotlib.backends.backend_pdf import PdfPages
from plot_param import parameters
parameters()

##############################################################################
################################ FUNCTIONS ###################################
def separate_sorted_df(df):
    '''
    sorts normalized df into combo types
    df: normalized, names not mapped
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

    bugs_R = df[df['Label_right'].str.contains('BUGS')]
    sep_dict['bug_abx'] = bugs_R[bugs_R['Label_left'].str.contains('abx')]

    media = df[df['Label_left'].str.contains('CAMHB')]
    sep_dict['media_pair'] = media[media['Label_right'].str.contains('CAMHB')]

    return sep_dict

##############################################################################
################################## INPUTS ####################################
norm_file = snakemake.input[0]
key_file = snakemake.input[1]

keys = pd.read_excel(key_file,index_col=0,sheet_name='ABX')
cpkeys = pd.read_excel(key_file,index_col=0,sheet_name='CP')
day = norm_file.split('/')[-1].split('_')[0]
chip = norm_file.split('/')[-1].split('_normalized')[0]
interim_folder = 'data/interim/'+day+'/'+chip
output_folder = 'data/output/'+day+'/'+chip
figure_folder='figures/'+day+'/'+chip

STRAIN = norm_file.split('/')[-1].split('_')[1]
CPP = norm_file.split('/')[-1].split('_')[2]

# create strain-abx dict from keys
ab1_key = dict(zip(keys.index, keys.ab1))
ab2_key = dict(zip(keys.index, keys.ab2))
kp1_key = dict(zip(keys.index, keys.kp1))
kp2_key = dict(zip(keys.index, keys.kp2))
pa1_key = dict(zip(keys.index, keys.pa1))
pa2_key = dict(zip(keys.index, keys.pa2))
strains = ['ab1','ab2','kp1','kp2','pa1','pa2']
strain_key = dict(zip(strains,[ab1_key,ab2_key,kp1_key,kp2_key,pa1_key,pa2_key]))
abxdict = strain_key[STRAIN]

# cp name dict
cppkeys = cpkeys[cpkeys.Blainey_cpp==CPP]
cppdict = dict(zip(cppkeys.Blainey_cp,cppkeys.Broad_Sample))

# save outputs
pp = PdfPages(figure_folder+'_stats_compiled.pdf')

##############################################################################
################################ PLOTTING ####################################
'''
TABLE: list of missing inputs
'''
missing_file = snakemake.input[2]
missing = pd.read_csv(missing_file,index_col=0)
if len(missing)>0:
    missing['abx_name'] = missing['missing_clusters'].map(abxdict)
    missing['cp_name'] = missing['missing_clusters'].map(cppdict)
    missing['chipID'] = chip
    missing.to_csv(missing_file)

    fig, ax =plt.subplots(1,1,figsize=(20,20))
    table = ax.table(cellText=missing.values,rowLabels=missing.index,colLabels=missing.columns,loc='top')
    ax.axis('tight')
    ax.axis('off')
    plt.savefig(pp, format='pdf',bbox_inches='tight',dpi=300)
else:
    pass

'''
SCATTER
    - Total # of droplets detected from droplets unprocessed (droplets_b)
    - Total # of droplets detected with 2 droplets (clustered)
    - Total # of microwells detected with 2 droplets (trimmed)
    - Total # of microwells detected with 2 droplets (qcfiltered)
    - Total # of microwells detected with 2 droplets (outlierfiltered)

'''
raw = pd.read_csv(snakemake.input[3],index_col=0)
clus = pd.read_csv(snakemake.input[4],index_col=0)
prefil = pd.read_csv(snakemake.input[5],index_col=0)
postfil = pd.read_csv(snakemake.input[6],index_col=0)
outlierfil = pd.read_csv(snakemake.input[7],index_col=0)

N_df = pd.DataFrame(data=[raw.shape[0],clus.shape[0],prefil.shape[0],postfil.shape[0],outlierfil.shape[0]],\
             index=['raw','clustered','prefilter','postfilter','outlier_removed'],columns=['numbers'])
N_df['percent'] = 0
N_df.loc[(N_df.index=='raw')|(N_df.index=='clustered'),'percent'] = N_df.numbers / N_df.loc['raw','numbers']
N_df.loc[(N_df.index=='prefilter')|(N_df.index=='postfilter')|(N_df.index=='outlier_removed'),\
         'percent'] = N_df.numbers / N_df.loc['prefilter','numbers']

fig,ax=plt.subplots(1,2,figsize=(20,12))
sns.barplot(data=N_df,x=N_df.index,y='numbers',ax=ax[0])
ax[0].set_title('Total No. of droplets or microwells')
sns.barplot(data=N_df,x=N_df.index,y='percent',ax=ax[1])
ax[1].set_title('Percent of droplets or microwells prefilter')
plt.suptitle(chip)
for ax in fig.axes:
    ax.tick_params(labelrotation=90)
plt.savefig(pp, format='pdf',bbox_inches='tight',dpi=300)

'''
SCATTER:
    - droplet counts, prefilter
    - droplet counts, postfilter
'''
def count_clusters(df):
    left_count = df['Label_left'].value_counts()
    right_count = df['Label_right'].value_counts()
    counts = left_count + right_count
    return counts

fig, ax = plt.subplots(3,1,figsize=(20,24))
ax[0].plot(count_clusters(prefil).values,'o')
ax[0].set_xticks(range(len(count_clusters(prefil).index.values)))
ax[0].set_xticklabels(count_clusters(prefil).index.values,size=10,rotation=90)
ax[0].set_ylim(0, 1.5*max(count_clusters(prefil).values))
ax[0].set_title(chip+ ' - cluster counts prefilter')
ax[1].plot(count_clusters(postfil).values,'o')
ax[1].set_xticks(range(len(count_clusters(postfil).index.values)))
ax[1].set_xticklabels(count_clusters(postfil).index.values,size=10,rotation=90)
ax[1].set_ylim(0, 1.5*max(count_clusters(postfil).values))
ax[1].set_title(chip+ ' - cluster counts postfilter')
ax[2].plot(count_clusters(outlierfil).values,'o')
ax[2].set_xticks(range(len(count_clusters(outlierfil).index.values)))
ax[2].set_xticklabels(count_clusters(outlierfil).index.values,size=10,rotation=90)
ax[2].set_ylim(0, 1.5*max(count_clusters(outlierfil).values))
ax[2].set_title(chip+ ' - cluster counts outliers removed')
for ax in fig.axes:
    ax.tick_params(labelrotation=90)
plt.savefig(pp, format='pdf',bbox_inches='tight',dpi=300)

'''
HISTOGRAM: control pair signal distributions
SCATTER:
    - Total # of replicates for CAMHB+CAMHB controls
    - Total # of replicates for CAMHB+BUG controls
    - Total # of replicates for BUG+BUG controls
'''
norm = pd.read_csv(norm_file,index_col=0)
bug_drops = norm[(norm.Label_left.str.contains('BUGS')) & \
(norm.Label_right.str.contains('BUGS'))]
bugmed_drops = norm[(norm.Label_left.str.contains('BUGS')) & \
(norm.Label_right.str.contains('CAMHB'))]
media_drops = norm[(norm.Label_left.str.contains('CAMHB')) & \
(norm.Label_right.str.contains('CAMHB'))]

fig,ax = plt.subplots(1,2,figsize=(20,6.3))
ax[0].hist(bug_drops.t0_norm2,label='BUGS',alpha=0.7)
ax[0].hist(bugmed_drops.t0_norm2,label='BUG+CAMHB',alpha=0.7)
ax[0].hist(media_drops.t0_norm2,label='CAMHB',alpha=0.7)
ax[0].set_title(chip+' - control distribution t0_norm2')
ax[0].set_xlabel('t0_norm2')
ax[0].legend(bbox_to_anchor=(1.1,1))

ctrl_N = [bug_drops.shape[0],bugmed_drops.shape[0],media_drops.shape[0]]
ax[1].bar(x=['BUGS','BUG+CAMHB','CAMHB'],height=ctrl_N)
ax[1].set_ylabel('No. of microwells')
ax[1].set_title(chip+' - control counts')
plt.savefig(pp, format='pdf',bbox_inches='tight',dpi=300)

'''
SCATTER: correlation plot
'''
norm_sort = norm.sort_values(by=['Label_left','Label_right'])
norm_even = norm_sort.iloc[::2,:].groupby(['Label_left','Label_right']).median()
norm_odd = norm_sort.iloc[1::2,:].groupby(['Label_left','Label_right']).median()
# drop combos that do not exist in both norm_even and norm_odd
norm_merged = norm_even.merge(norm_odd,how='outer',left_index=True,right_index=True).dropna()

corr, _ = pearsonr(norm_merged.t0_norm2_x, norm_merged.t0_norm2_y)

plt.figure(figsize=(20,20))
sns.scatterplot(data=norm_merged,x='t0_norm2_x',y='t0_norm2_y')
plt.title(chip+f' - correlation t0 norm2, R = {corr:.4f}')
plt.savefig(pp, format='pdf',bbox_inches='tight',dpi=300)

'''
BARPLOT: standard zprimes
'''
def z_factor(pos_,neg_,pos_std,neg_std):
    z = 1 - (3*(pos_std + neg_std) / (abs(pos_ - neg_)))
    return z

def calc_zfactor(df,negctrl_l,negctrl_r,posctrl_l,posctrl_r):
    '''
    negctrl_: 'CAMHB'
    posctrl_: 'BUGS'
    '''
    neg_SE = df[(df.Label_left.str.contains(negctrl_l)) & \
    (df.Label_right.str.contains(negctrl_r))]['t0_norm2_SE'].mean()
    neg_med = df[(df.Label_left.str.contains(negctrl_l)) & \
    (df.Label_right.str.contains(negctrl_r))]['t0_norm2'].median()

    pos_SE = df[(df.Label_left.str.contains(posctrl_l)) & \
    (df.Label_right.str.contains(posctrl_r))]['t0_norm2_SE'].mean()
    pos_med = df[(df.Label_left.str.contains(posctrl_l)) & \
    (df.Label_right.str.contains(posctrl_r))]['t0_norm2'].median()
    z_ = z_factor(pos_med,neg_med,pos_SE,neg_SE)
    return z_

bliss_dict = separate_sorted_df(norm)
lowest = bliss_dict['abx_media'].groupby('Label_left').median().reset_index().sort_values(by='norm_growth').Label_left[0]
z_2X = calc_zfactor(norm,'CAMHB','CAMHB','BUGS','BUGS')
z_1X = calc_zfactor(norm,'CAMHB','CAMHB','BUGS','CAMHB')

z_low_2X = calc_zfactor(norm,lowest,'CAMHB','BUGS','BUGS')
z_low_1X = calc_zfactor(norm,lowest,'CAMHB','BUGS','CAMHB')

plt.figure(figsize=(20,20))
zprimes_labels = ['2XBUG','1XBUG','2XBUG_'+lowest,'1XBUG_'+lowest]
zprimes = [z_2X,z_1X,z_low_2X,z_low_1X]
plt.bar(x=zprimes_labels,height=zprimes)
plt.axhline(y=0.4, color='k', linestyle=':')
plt.xticks(rotation=90)
plt.ylim([0,1.2])
plt.title(chip+' - standard z primes')
plt.savefig(pp, format='pdf',bbox_inches='tight',dpi=300)
zprime_st = pd.DataFrame(np.asarray([zprimes_labels,zprimes]).T,columns=['zprime','value'])
zprime_st['chipID'] = chip
zprime_st.to_csv(snakemake.output[0])

'''
BARPLOT: abx zprimes
'''
def abx_zfactor(abx_X,bug_X):
    abx_df = bliss_dict[abx_X].groupby('Label_left').median().reset_index()
    abx_df['abx'] = abx_df['Label_left'].apply(lambda x: x[:-1])
    abx_df['conc'] = abx_df['Label_left'].apply(lambda x: x[-1])
    abx_a = abx_df[abx_df.conc=='a']
    a_med_dict = dict(zip(abx_a.abx,abx_a.t0_norm2))
    a_SE_dict = dict(zip(abx_a.abx,abx_a.t0_norm2_SE))
    abx_a['a_t0norm2'] = abx_a['abx'].map(a_med_dict)
    abx_a['a_t0norm2_SE'] = abx_a['abx'].map(a_SE_dict)
    abx_a['bug_med'] = bliss_dict[bug_X].median()['t0_norm2']
    abx_a['bug_SE'] = bliss_dict[bug_X].median()['t0_norm2_SE']
    abx_a['abx_zfactor'] = 1 - (3*(abx_a['bug_SE'] + abx_a['a_t0norm2_SE']) / \
                                      (abs(abx_a['bug_med'] - abx_a['a_t0norm2'])))
    abx_a.replace({'abx':abxdict},inplace=True)
    return abx_a

bug1x_abx = abx_zfactor('abx_media','bug_media')
bug2x_abx = abx_zfactor('abx_bugs','bug_pair')

fig,ax = plt.subplots(2,1,figsize=(20,20))
sns.barplot(x='abx', y='abx_zfactor', data=bug1x_abx, ax=ax[0])
ax[0].set_title(chip+' - abx specific zprime 1XBUG')
ax[0].set_ylim([0,1.2])
ax[0].axhline(y=0.4, color='k', linestyle=':')
sns.barplot(x='abx', y='abx_zfactor', data=bug2x_abx, ax=ax[1])
ax[1].set_title(chip+' - abx specific zprime 2XBUG')
ax[1].set_ylim([0,1.2])
ax[1].axhline(y=0.4, color='k', linestyle=':')
for ax in fig.axes:
    ax.tick_params(labelrotation=90)
plt.savefig(pp, format='pdf',bbox_inches='tight',dpi=300)
zprime_abx = bug1x_abx.merge(bug2x_abx,how='outer',on='Label_left',suffixes=['_1X','_2X'])
zprime_abx['chipID'] = chip
zprime_abx.to_csv(snakemake.output[1])

'''
HISTOGRAM: replicate numbers for abx-abx, abx-camhb, cp-bugs, and abx-cp combos
'''
abx_abx = bliss_dict['abx_pair'].groupby(['Label_left','Label_right']).count()['t0_norm2']
abx_med = bliss_dict['abx_media'].groupby(['Label_left','Label_right']).count()['t0_norm2']
bug_cp = bliss_dict['bug_cp'].groupby(['Label_left','Label_right']).count()['t0_norm2']
abx_cp = bliss_dict['abx_cp'].groupby(['Label_left','Label_right']).count()['t0_norm2']

plt.figure(figsize=(20,10))
plt.hist(abx_abx,density=True,alpha=0.5,label='abx-abx')
plt.hist(abx_med,density=True,alpha=0.5,label='abx-media')
plt.hist(bug_cp,density=True,alpha=0.5,label='bug-cp')
plt.hist(abx_cp,density=True,alpha=0.5,label='abx-cp')
plt.xlabel('No. of microwells')
plt.title(chip+' - No. of microwells')
plt.legend(title='combo type',bbox_to_anchor=(1.1,1))
plt.savefig(pp, format='pdf',bbox_inches='tight',dpi=300)

'''
SCATTER: droplet positions on chip
'''
pre_post = pd.read_csv(snakemake.input[8],index_col=0)

fig, axes = plt.subplots(2,1,figsize=(18,36))
axes[0].plot(pre_post['Pre_GlobalX'],-pre_post['Pre_GlobalY'],'.',ms=1)
axes[0].set_title(chip+' - wells detected')

axes[1].plot(pre_post['Post_GlobalX'],-pre_post['Post_GlobalY'],'.',ms=1)
axes[1].set_title(chip+' - wells detected - edge wells removed')
plt.tight_layout()
plt.savefig(pp, format='pdf',bbox_inches='tight',dpi=300)

'''
HISTOGRAM: bug-cp growth (inhibition)
'''
cp_only = pd.read_csv(snakemake.input[9],index_col=0)

f, (ax_hist, ax_box) = plt.subplots(2, sharex=True, figsize=(20,10),\
gridspec_kw={"height_ratios": (.85, .15)})
sns.distplot(cp_only.norm_growth, ax=ax_hist)
sns.boxplot(cp_only.norm_growth, ax=ax_box)
ax_box.set(yticks=[])
sns.despine(ax=ax_hist)
sns.despine(ax=ax_box, left=True)
ax_box.set_xlim([-1,2])
ax_hist.set_xlim([-1,2])
plt.suptitle(chip+' - compound only activity, normalized growth')
plt.savefig(pp, format='pdf',bbox_inches='tight',dpi=300)

f, (ax_hist, ax_box) = plt.subplots(2, sharex=True, figsize=(20,10),\
gridspec_kw={"height_ratios": (.85, .15)})
sns.distplot(cp_only.norm_growth_inh, ax=ax_hist)
sns.boxplot(cp_only.norm_growth_inh, ax=ax_box)
ax_box.set(yticks=[])
sns.despine(ax=ax_hist)
sns.despine(ax=ax_box, left=True)
ax_box.set_xlim([-1,2])
ax_hist.set_xlim([-1,2])
plt.suptitle(chip+' - compound only activity, normalized growth inhibition')
plt.savefig(pp, format='pdf',bbox_inches='tight',dpi=300)

'''
ERRORBAR: abx curves on single axis
'''
abx_only = pd.read_csv(snakemake.input[10],index_col=0)

# num = abx_only.shape[0]
# # clr_map = plt.cm.viridis # LinearSegmentedColormap
# # N_clrs = min(clr_map.N,num)
# # map_clrs = [clr_map(int(x*clr_map.N/N_clrs)) for x in range(N_clrs)]
# plt.figure(figsize=(20,3))
# for i in np.arange(0, num, 3):
#     plt.errorbar(x=abx_only.abx_name_conc[i:i+3], y=abx_only.norm_growth[i:i+3], \
# 	yerr=abx_only.norm_growth_SE[i:i+3]) #, c=map_clrs[i])
#
# unique = np.unique(abx_only.abx_name)
# plt.xticks(np.arange(0, num, 3), unique, rotation=45)
# plt.ylabel('median relative growth')
# plt.xlabel('antibiotics in descending concentration order')
# plt.title(chip+' - all abx curves')
# plot all abx on single axeis
num = abx_only.abx_name.unique().shape[0]
# clr_map = plt.cm.viridis # LinearSegmentedColormap
# N_clrs = min(clr_map.N,num)
# map_clrs = [clr_map(int(x*clr_map.N/N_clrs)) for x in range(N_clrs)]
plt.figure(figsize=(20,3))
# for i in np.arange(0, num, 3):
#     plt.errorbar(x=abxMED.abx_name_conc[i:i+3], y=abxMED.norm_growth[i:i+3], \
# 	yerr=abxMED.norm_growth_SE[i:i+3], c=map_clrs[i])
# changed
xticks_ = [0]
val = 0
for i in enumerate(abx_only.abx_name.unique()):
    plt.errorbar(x=abx_only.loc[(abx_only.abx_name==i[1]),'abx_name_conc'], y=abx_only.loc[(abx_only.abx_name==i[1]),'norm_growth'], \
	yerr=abx_only.loc[(abx_only.abx_name==i[1]),'norm_growth_SE'])#, c=map_clrs[i[0]])
    val+=abx_only[abx_only.abx_name==i[1]].shape[0]
    xticks_.append(val)

unique = np.unique(abx_only.abx_name)
# plt.xticks(np.arange(0, num, 3), unique, rotation=45)
plt.xticks(xticks_[:-1], unique, rotation=45)
plt.ylabel('median relative growth')
plt.xlabel('antibiotics in descending concentration order')
plt.title(chip+' - all abx curves')
plt.savefig(pp, format='pdf',bbox_inches='tight',dpi=300)

'''
HEATMAP: abx-cp summed bliss score + se
'''
abx_abx = pd.read_csv(snakemake.input[11],index_col=0)

def matrixize(df,ind,col,val):
    '''
    Makes melted df into a filled square matrix
    '''
    df_pivot = df.pivot(index=ind,columns=col, values=val)
    df_pivot = df_pivot.fillna(0)
    df_fill = df_pivot.values + df_pivot.values.T
    np.fill_diagonal(df_fill,0.5*np.diag(df_fill))
    df_mat = pd.DataFrame(data=df_fill,index=df_pivot.index,columns=df_pivot.columns)
    return df_mat

def matrix_mask(df):
    '''
    Creates mask for heatmap of square matrix
    '''
    mask = np.zeros_like(df, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True
    return mask

# matrixizing data
df_bliss = matrixize(abx_abx,ind='abx_name1',col='abx_name2',val='bliss_med')
mask_bliss = matrix_mask(df_bliss)
df_se = matrixize(abx_abx,ind='abx_name1',col='abx_name2',val='bliss_se')
mask_se = matrix_mask(df_se)

max = np.quantile(abs(abx_abx.bliss_med),.98)
fig, ax = plt.subplots(1,2,figsize=(20,8))
ax[0].set_title(chip+' - abx-abx summed bliss scores')
ax[1].set_title(chip+' - abx-abx summed bliss SE')
sns.heatmap(df_bliss,cmap='coolwarm',mask=mask_bliss,vmin=-max,vmax=max,\
cbar_kws={'pad':0.01,'aspect':50},ax=ax[0])
sns.heatmap(df_se,cmap='coolwarm',mask=mask_se,cbar_kws={'pad':0.01,'aspect':50},ax=ax[1])
plt.savefig(pp, format='pdf',bbox_inches='tight',dpi=300)

'''
HEATMAP: abx-cp norm growth + count
'''
# extract media and compound - abx combos only
abx_medcp = pd.concat([bliss_dict['abx_media'],bliss_dict['abx_cp']])

# mapping
abx_medcp['abx'] = abx_medcp.Label_left.apply(lambda x: x[:-1])
abx_medcp['abx_name'] = abx_medcp['abx'].map(abxdict)
abx_medcp['abx_conc'] = abx_medcp.Label_left.apply(lambda x: x[-1])
abx_medcp['abx_name_conc'] = abx_medcp.abx_name+'_'+abx_medcp.abx_conc
abx_medcp.replace({'Label_right': cppdict},inplace=True)

# median
abx_medcp_ = abx_medcp.groupby(['abx_name_conc',\
'Label_right']).median().reset_index().sort_values(by='abx_name_conc')
abx_medcp_mat = abx_medcp_.pivot(index='abx_name_conc',columns='Label_right', values='norm_growth')

# counts
abx_medcp_ct = abx_medcp.groupby(['abx_name_conc',\
'Label_right']).count().reset_index().sort_values(by='abx_name_conc')
abx_medcp_ct_ = abx_medcp_ct.pivot(index='abx_name_conc',columns='Label_right', values='norm_growth')

# plotting abx-cp growth norm heatmap
fig, ax = plt.subplots(1,2,figsize=(20,8))
ax[0].set_title(chip+' - abx-cp normalized growth')
ax[1].set_title(chip+' - abx-cp counts')
sns.heatmap(abx_medcp_mat,vmin=0,vmax=1.5,cmap='coolwarm',cbar_kws={'pad':0.01,'aspect':50},ax=ax[0])
sns.heatmap(abx_medcp_ct_,cmap='coolwarm',annot=True,cbar_kws={'pad':0.01,'aspect':50},annot_kws={"size":5},ax=ax[1])
plt.savefig(pp, format='pdf',bbox_inches='tight',dpi=300)

# close output
pp.close()

##############################################################################
################################ CONCAT PDF ##################################
# _compiled_chip_hitstats
pdfs = [snakemake.input[12],figure_folder+'_stats_compiled.pdf']

merger = PdfFileMerger()

for pdf in pdfs:
    merger.append(pdf)

merger.write(figure_folder+'_compiled_chip_hitstats.pdf')
merger.close()

# _compiled_all_figures
image_files = glob.glob(figure_folder+'*.png')
imagelist = []
for i in image_files:
    image1 = Image.open(i)
    im1 = image1.convert('RGB')
    imagelist.append(im1)

imagelist[0].save(snakemake.output[2],save_all=True,\
append_images=imagelist)
