'''
Script: abx_cp.py
Last update: 2021 March, mzhu

Description:
    1. abx-cp bliss score and se
    2. abx-cp null bliss score and se
    3. t-test for signficance
    4. chip hits based on synergy ctrl bliss and se (pval calc with same abx-cp null)
    5. plot chip hits info
Inputs:
    1. data/output/{DAY}/{chip}_normalized.csv
    2. data/output/{DAY}/{chip}_abxabx_synctrl.xlsx
    3. Keys-*.xlsx
    4. qc_params.csv
Outputs (csv):
    1. data/output/{DAY}/{chip}_abxcp_droplet_data.csv
    2. data/output/{DAY}/{chip}_abxcp_negctrl_droplet_data.csv
    3. data/output/{DAY}/{chip}_abxcp_negctrl_summary_data.csv
    4. data/output/{DAY}/{chip}_abxcp_negctrl_summed_summary_data.csv
    5. data/output/{DAY}/{chip}_abxcp_summed_summary_data.csv
    6. data/output/{DAY}/{chip}_abxcp_summary_data.csv
    7. data/output/{DAY}/{chip}_abxcp_chip_hit.csv
Outputs (figures):
    1. figures/{DAY}/{chip}_abxcp_normgrowth_heatmap.png
    2. figures/{DAY}/{chip}_abxcp_bliss_distribution.png
    3. figures/{DAY}/{chip}_abxcp_null_distribution.png
    4. figures/{DAY}/{chip}_abxcp_bliss_volcano.png
    5. figures/{DAY}/{chip}_abxcp_chiphits.pdf
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
import math as mt
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import multiprocessing as mp
from pathos.multiprocessing import ProcessingPool as Pool
import os, time
import warnings
warnings.filterwarnings('ignore')
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.image as image
from matplotlib.backends.backend_pdf import PdfPages
from plot_param import parameters
parameters()

##############################################################################
################################# FUNCTIONS ##################################
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

def bootstrap_median(df):
    ''' resample from all dataframe rows '''
    return np.median(np.random.choice(df, len(df)))

def bliss_calc(abx_cp, ctrl_a, ctrl_c, parallel=False, save=False):
    '''
    Calculates Bliss scores for given abx-cp pair.

    Inputs:
    :abx_cp: (DataFrame) abx-cp combo data, normalized to bug-media
    :ctrl_a: (DataFrame) abx-camhb combo data, normalized to bug-media
    :ctrl_c: (DataFrame) cp-bug combo data, normalized to bug-media

    Outputs:
    if parallel: just the bliss scores
    if non-parallel: entire DataFrame with E_a, E_c, and bliss columns added
    '''

    med_agg = bootstrap_median if parallel else 'median'

    abg = ctrl_a[['Label_left','norm_growth_inh']].groupby('Label_left')
    control_a = abg.aggregate(med_agg)

    cbg = ctrl_c[['Label_right','norm_growth_inh']].groupby('Label_right')
    control_c = cbg.aggregate(med_agg)

    data_df = abx_cp[['Label_left', 'Label_right', 'norm_growth_inh']].groupby(['Label_left', 'Label_right'])
    data_df = data_df.aggregate(med_agg)

    synergy_df = data_df.merge(control_a,left_index=True,right_index=True)
    synergy_df.rename(columns={'norm_growth_inh_y':'E_a'}, inplace=True)

    synergy_df = synergy_df.merge(control_c,left_index=True,right_index=True)
    synergy_df.rename(columns={'norm_growth_inh':'E_c', 'norm_growth_inh_x': 'E_ac'}, inplace=True)
    bliss = synergy_df.copy()
    bliss['bliss'] = (bliss.E_ac - (bliss.E_a + bliss.E_c - bliss.E_a*bliss.E_c)).astype('float')

    if parallel:
    	return bliss.loc[:, 'bliss'], bliss.loc[:,'E_a'],  bliss.loc[:,'E_c'],  bliss.loc[:,'E_ac']
    else:
    	return bliss

def bliss_parallel(arglist):
    return bliss_calc(*arglist)

##############################################################################
norm_file = snakemake.input[0]
syn_ctrl = pd.read_csv(snakemake.input[1],index_col=0)
key_file = snakemake.input[2]
keys = pd.read_excel(key_file,index_col=0,sheet_name='ABX')
cpkeys = pd.read_excel(key_file,index_col=0,sheet_name='CP')

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
CPP = norm_file.split('/')[-1].split('_')[2]

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

abxdict = strain_key[STRAIN]

# CREATING CP DICTIONARY FROM KEYS
cppkeys = cpkeys[cpkeys.Blainey_cpp==CPP]
cppdict = dict(zip(cppkeys.Blainey_cp,cppkeys.Broad_Sample))

##############################################################################
######################## ABX CP HEATMAP ######################################
# create dictionary of df combo types
bliss_dict = separate_sorted_df(norm)

# extract media and compound - abx combos only
abx_medcp = pd.concat([bliss_dict['abx_media'],bliss_dict['abx_cp']])

# map abx and cp names
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

# plotting heatmaps
fig, ax = plt.subplots(1,2,figsize=(20,8))
ax[0].set_title(chip+' - abx-cp normalized growth')
ax[1].set_title(chip+' - abx-cp counts')
sns.heatmap(abx_medcp_mat,cmap='coolwarm',cbar_kws={'pad':0.01,'aspect':50},ax=ax[0])
sns.heatmap(abx_medcp_ct_,cmap='coolwarm',annot=True,cbar_kws={'pad':0.01,\
'aspect':50},annot_kws={"size":5},ax=ax[1])
plt.savefig(figure_folder+'_abxcp_normgrowth_heatmap.png')

##############################################################################
######################## ABX CP BLISS ########################################
# parallelized bliss scoring
# create the multiprocessing pool
cores = mp.cpu_count()//int(POOL_NUM)
pool = Pool(cores)

# process the dataFrame by mapping function to each df across the pool
bliss_pool, bliss_Ea, bliss_Ec, bliss_Eac = zip(*pool.map(bliss_parallel,
                    ([bliss_dict['abx_cp'],
                    bliss_dict['abx_media'],
                    bliss_dict['bug_cp'],
                    i] for i in range(1,1001))))

# close down the pool and join
pool.close()
pool.join()
pool.clear()

# calc bliss and se for abx-cp
bliss_boot = pd.concat(bliss_pool, axis=1)
bliss_boot_Ea = pd.concat(bliss_Ea, axis=1)
bliss_boot_Ec = pd.concat(bliss_Ec, axis=1)
bliss_boot_Eac = pd.concat(bliss_Eac, axis=1)

def boot_med_SE(df_boot):
    '''
    Inputs:
    :df_boot: bootstrapped df

    Outputs:
    :_med_dict: comboID --> median vals dict
    :_se_dict: comboID --> std vals dict
    '''
    _med_bliss = df_boot.transpose().median().reset_index()
    _med_bliss['comboID'] = list(zip(_med_bliss.Label_left,_med_bliss.Label_right))
    _med_dict = dict(zip(_med_bliss.comboID,_med_bliss.iloc[:,2]))

    _se_bliss = df_boot.transpose().std().reset_index()
    _se_bliss['comboID'] = list(zip(_se_bliss.Label_left,_se_bliss.Label_right))
    _se_dict = dict(zip(_se_bliss.comboID,_se_bliss.iloc[:,2]))
    return _med_dict, _se_dict

# get bootstrapped value dict
bliss_med_dict, bliss_se_dict = boot_med_SE(bliss_boot)
Ea_med_dict, Ea_se_dict = boot_med_SE(bliss_boot_Ea)
Ec_med_dict, Ec_se_dict = boot_med_SE(bliss_boot_Ec)
Eac_med_dict, Eac_se_dict = boot_med_SE(bliss_boot_Eac)

# add bliss and se to abx-cp df
abxcp = bliss_dict['abx_cp']
abxcp['comboID'] = list(zip(abxcp.Label_left,abxcp.Label_right))
abxcp['bliss_med'] = abxcp['comboID'].map(bliss_med_dict)
abxcp['bliss_se'] = abxcp['comboID'].map(bliss_se_dict)
abxcp['Ea_med'] = abxcp['comboID'].map(Ea_med_dict)
abxcp['Ea_se'] = abxcp['comboID'].map(Ea_se_dict)
abxcp['Ec_med'] = abxcp['comboID'].map(Ec_med_dict)
abxcp['Ec_se'] = abxcp['comboID'].map(Ec_se_dict)
abxcp['Eac_med'] = abxcp['comboID'].map(Eac_med_dict)
abxcp['Eac_se'] = abxcp['comboID'].map(Eac_se_dict)

# map abx and cp names
abxcp['abx'] = abxcp.Label_left.apply(lambda x: x[:-1])
abxcp['abx_name'] = abxcp['abx'].map(abxdict)
abxcp['cp_name'] = abxcp['Label_right'].map(cppdict)
abxcp['abx_conc'] = abxcp.Label_left.apply(lambda x: x[-1])
abxcp['abx_name_conc'] = abxcp.abx_name+'_'+abxcp.abx_conc

# calc tstat
abxcp['bliss_tstat'] = abxcp.bliss_med / abxcp.bliss_se

# save droplet data
abxcp['chipID'] = chip
abxcp.to_csv(snakemake.output[0])

# calc and save summary data
abxcp_med = abxcp[['Label_right','abx_name','cp_name','abx_name_conc',\
'bliss_med','bliss_se','Ea_med','Ea_se','Ec_med','Ec_se','Eac_med','Eac_se',\
'bliss_tstat']].groupby(['Label_right','abx_name','cp_name','abx_name_conc'\
]).median().reset_index()
abxcp_med['chipID'] = chip

# calc summed bliss scores
abxcp_sum = abxcp_med[['Label_right','abx_name','cp_name','bliss_med',\
'bliss_se']].groupby(['Label_right','abx_name','cp_name']).sum().reset_index()
abxcp_sum['bliss_tstat'] = abxcp_sum.bliss_med / abxcp_sum.bliss_se
abxcp_sum['chipID'] = chip

# plotting summed bliss score distribution
plt.figure()
plt.title(chip+' - abx-cp distribution')
sns.distplot(abxcp_sum['bliss_med'])
plt.savefig(figure_folder+'_abxcp_bliss_distribution.png')

##############################################################################
######################## NEGATIVE CONTROLS ###################################
def ctrl_bliss_calc(abxmed1, abxmed2, bugmed, parallel=False, save=False):
    '''
    Calculates Bliss scores for abx cp null distribution.

    Inputs:
    :abxmed1: (DataFrame) abx-media combo data, normalized to bug-media
    :abxmed2: (DataFrame) abx-media combo data, normalized to bug-media
    :bugmed: (DataFrame) bug-media combo data, normalized to bug-media

    Outputs:
    if parallel: just the bliss scores
    if non-parallel: entire DataFrame with E_a, E_b, and bliss columns added
    '''

    med_agg = bootstrap_median if parallel else 'median'

    abg = abxmed1[['Label_left','norm_growth_inh']].groupby('Label_left')
    control_a = abg.aggregate(med_agg) # controls grouped by abx

    cbg = bugmed.norm_growth_inh.values
    control_b = np.median(np.random.choice(cbg,len(cbg)))

    data_df = abxmed2[['Label_left','norm_growth_inh']].groupby(['Label_left'])
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

# extract negative controls
neg_Eab = bliss_dict['abx_media'].sort_values(by=['Label_left','Label_right']).iloc[::2,:]
neg_Ea = bliss_dict['abx_media'].sort_values(by=['Label_left','Label_right']).iloc[1::2,:]
neg_Eb = bliss_dict['bug_media']

# parallelized bliss scoring
# create the multiprocessing pool
cores = mp.cpu_count()//int(POOL_NUM)
pool = Pool(cores)

# process the dataFrame by mapping function to each df across the pool
ctrl_bliss_pool = pool.map(ctrl_bliss_parallel,([neg_Eab,neg_Ea,neg_Eb,i] for i in range(1,1001)))

# close down the pool and join
pool.close()
pool.join()
pool.clear()

# calc median bliss and se for controls
ctrl_bliss_boot = pd.concat(ctrl_bliss_pool, axis=1)
ctrl_med_bliss = ctrl_bliss_boot.transpose().median().reset_index()
ctrl_med_dict = dict(zip(ctrl_med_bliss.Label_left,ctrl_med_bliss.iloc[:,1]))

ctrl_se_bliss = ctrl_bliss_boot.transpose().std().reset_index()
ctrl_se_dict = dict(zip(ctrl_se_bliss.Label_left,ctrl_se_bliss.iloc[:,1]))

# add median bliss and se to abx-media df
abxmed = bliss_dict['abx_media']
abxmed['bliss_med'] = abxmed['Label_left'].map(ctrl_med_dict)
abxmed['bliss_se'] = abxmed['Label_left'].map(ctrl_se_dict)

# map abx names
abxdict = strain_key[STRAIN]
abxmed['abx'] = abxmed.Label_left.apply(lambda x: x[:-1])
abxmed['abx_name'] = abxmed['abx'].map(abxdict)
abxmed['abx_conc'] = abxmed.Label_left.apply(lambda x: x[-1])
abxmed['abx_name_conc'] = abxmed.abx_name+'_'+abxmed.abx_conc

# calc tstat
abxmed['bliss_tstat'] = abxmed.bliss_med / abxmed.bliss_se

# calc median bliss score and se
abxmed_med = abxmed[['abx_name_conc','bliss_med',\
'bliss_se','bliss_tstat']].groupby('abx_name_conc').median().reset_index()
abxmed_med['abx_name'] = abxmed_med.abx_name_conc.apply(lambda x: x.split('_')[0])

# calc 3 sum bliss scores and se
abxmed_sum = abxmed_med[['abx_name','bliss_med',\
'bliss_se']].groupby('abx_name').sum()

# calc summed tstat
abxmed_sum['bliss_tstat'] = abxmed_sum.bliss_med / abxmed_sum.bliss_se

# add chip info
abxmed['chipID'] = chip
abxmed_med['chipID'] = chip
abxmed_sum['chipID'] = chip

# save outputs
abxmed.to_csv(snakemake.output[1])
abxmed_med.to_csv(snakemake.output[2])
abxmed_sum.to_csv(snakemake.output[3])

##############################################################################
############################# T-TEST ########################################
# fit null to t-distribution
null = abxmed_sum['bliss_tstat'].fillna(0).values.tolist()
tdof, tloc, tscale = stats.t.fit(null, floc = 0.)
null_fit = stats.t.rvs(tdof, loc=tloc, scale=tscale, size=10000)
rv = stats.t(tdof, loc=tloc, scale=tscale)

# plotting null tstat score distribution
a = np.arange(-5,5,0.1)
plt.figure()
plt.title(chip+' - abx-cp null distribution')
plt.plot(a, rv.pdf(a), '-', label='fitted data')
plt.hist(abxmed_sum['bliss_tstat'],density=True,label='empirical data')
plt.legend(bbox_to_anchor=(1.25,1))
plt.savefig(figure_folder+'_abxcp_null_distribution.png')

def t_test(tstat):
    ''' Get p value from t test '''
    global tdof, tloc, tscale
    return stats.t.sf(abs(tstat), tdof, loc=tloc, scale=tscale)

# add in synergy control
abxcp_sum['synctrl'] = 'cp'
syn_ctrl['synctrl'] = syn_ctrl.abx1_name+syn_ctrl.abx2_name
syn_ctrl.rename(columns={'abx1_name':'abx_name','abx2_name':'cp_name'},inplace=True)
abxcp_sum = pd.concat([syn_ctrl[['bliss_med','bliss_se','abx_name','cp_name',\
'bliss_tstat','chipID','synctrl']],abxcp_sum])

# calc p value from t test
abxcp_sum['bliss_pval'] = abxcp_sum['bliss_tstat'].apply(t_test)
abxcp_sum['bliss_log10pval'] = -np.log10(abxcp_sum['bliss_pval'])
abxcp_sum.sort_values(by='bliss_med',ascending=False,inplace=True)
abxcp_sum['chip_hit'] = 0
abxcp_sum.loc[(abxcp_sum.bliss_med>=syn_ctrl.bliss_med.values[0]),'chip_hit'] = 1

# save output
abxcp_sum.to_csv(snakemake.output[4])

# plot volcano
plt.figure()
sns.scatterplot(data=abxcp_sum,x='bliss_med',y='bliss_log10pval',hue='synctrl',alpha=0.7)
plt.title(chip+' - abx-cp volcano plot')
plt.legend(bbox_to_anchor=(1.25,1))
plt.savefig(figure_folder+'_abxcp_bliss_volcano.png')

##############################################################################
############################# CHIP HITS #######################################
# chip hit plotting
chiphit = abxcp_sum[abxcp_sum.chip_hit==1]
chemprop = pd.read_excel('crl_props.xlsx').fillna('EMPTY')

# creating combo and cp hit dict
chiphit['combohitID'] = chiphit.abx_name+chiphit.Label_right
hit_dict = dict(zip(chiphit['combohitID'],chiphit['chip_hit']))
cphit_dict = dict(zip(chiphit['cp_name'],chiphit['chip_hit']))

# mapping combo and cp hits onto abxcp_med
abxcp_med['comboID'] = abxcp_med.abx_name+abxcp_med.Label_right
abxcp_med['chip_cphit'] = 0
abxcp_med['chip_cphit'] = abxcp_med['cp_name'].map(cphit_dict)
abxcp_med['chip_combohit'] = 0
abxcp_med['chip_combohit'] = abxcp_med['comboID'].map(hit_dict)
abxcp_med_cphit = abxcp_med[abxcp_med.chip_cphit==1]

# plotting chip hits
if abxcp_med_cphit.shape[0]>0:
    pp = PdfPages(snakemake.output[7])
    for cp in abxcp_med_cphit.cp_name.unique():
        try:
            cp_df = abxcp_med_cphit[abxcp_med_cphit.cp_name==cp]
            chemprop_cp = chemprop[chemprop.ID==cp][['ID','MOLECULAR_WEIGHT','ALOGP','SMILES']].T
            chemprop_cp.columns = chemprop_cp.iloc[0]
            chemprop_cp = chemprop_cp.iloc[1:,:]
            # chem prop table and structure
            fig_t, ax_t =plt.subplots(1,figsize=(10,10),sharex=True)
            table = ax_t.table(cellText=chemprop_cp.values,rowLabels=chemprop_cp.index,\
            colLabels=chemprop_cp.columns,loc='left')
            smiles = chemprop_cp.loc['SMILES',:][0]
            m = Chem.MolFromSmiles(smiles)
            im = Draw.MolToImage(m,size=(500,500))
            ax_t.imshow(im,aspect='equal')
            ax_t.set_xticks([])
            ax_t.set_yticks([])
            plt.savefig(pp, format='pdf',bbox_inches='tight',dpi=300)
        except:
            pass

        # plotting Ea, Ec, Eac values
        abxs = cp_df.abx_name.unique()
        figshape = mt.ceil(np.sqrt(len(abxs)))
        fig,ax = plt.subplots(figshape,figshape,figsize=(5*figshape,5*figshape))
        for abx,ax in list(zip(abxs,ax.ravel())):
            abx_cp_df = cp_df[cp_df.abx_name==abx].sort_values(by='abx_name_conc')
            abxhit_dict = dict(zip(abx_cp_df.abx_name,abx_cp_df.chip_combohit))
            ax.errorbar(x=abx_cp_df['abx_name_conc'],y=abx_cp_df['Eac_med'], \
                        yerr=abx_cp_df['Eac_se'], label='abx-cp')
            ax.errorbar(x=abx_cp_df['abx_name_conc'],y=abx_cp_df['Ea_med'], \
                        yerr=abx_cp_df['Ea_se'], label='abx-media')
            ax.errorbar(x=abx_cp_df['abx_name_conc'],y=abx_cp_df['Ec_med'], \
                        yerr=abx_cp_df['Ec_se'], marker='x',label='bug-cp')
            ax.legend()
            ax.set_ylim([-0.1,2])
            ax.set_title(abx+' Hit: '+str(abxhit_dict[abx]))
        plt.suptitle(cp+' - Normalized growth inhibition',y=0.9)
        plt.savefig(pp, format='pdf',bbox_inches='tight',dpi=300)

    # close output
    pp.close()
else: # generates a blank pdf if no chip hits
    pp = PdfPages(snakemake.output[7])
    plt.figure()
    plt.savefig(pp, format='pdf',bbox_inches='tight',dpi=300)
    pp.close()

# save outputs
abxcp_med.to_csv(snakemake.output[5])
chiphit.to_csv(snakemake.output[6])
