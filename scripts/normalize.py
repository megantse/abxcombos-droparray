'''
Script: normalize.py
Last update: 2024 Feb, mtse

Description: NORMALIZE t0_norm2 to 1X or 2X BUGS
Input:
    1. data/output/{DAY}/{chip}_qcfiltered.csv # if OUTLIER==0
    2. data/output/{DAY}/{chip}_outlierfiltered.csv # if OUTLIER==1
    3. qc_params.csv
Output:
    1. data/output/{DAY}/{chip}_normalized.csv
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
################################# FUNCTIONS ##################################
def normalize(df, val=None, tp='t0'):
	'''
    Inputs:
	:df: df you to normalize
	:val: number if not normalizing
    :tp: timepoint name

    Outputs:
    :_df: normalized df
	'''
	_df = df.copy()
	if val==0:
		_df['norm_growth'] = 0
		_df['norm_growth_inh'] = 0
	elif val==1:
		_df['norm_growth'] = 1
		_df['norm_growth_inh'] = 0
	else:
		val_list = [val]*len(_df)
		_df['norm_growth'] = _df[tp+'_norm2'] / val_list
		_df['norm_growth_inh'] = (val_list - _df[tp+'_norm2']) / val_list
		_df.loc[(_df['norm_growth_inh'] <= 0), 'norm_growth_inh'] = 0

	return _df

def get_growth_change(df, tp='t0'):
    '''
    Calculates normalized growth change and % growth inhibition.
    Normalization is against growth-changed

    Inputs:
    :df: (DataFrame) droplets data, normalized against media-only drops.
    Labels must be derep sorted.
    :tp: (list) List of timepoint labels

    Outputs: Droplets DataFrame with normalized values
    '''

    _df = df.copy()

    # split
    _abx = _df[_df['Label_left'].str.contains('abx')]
    _med = _df[_df['Label_left'].str.contains('CAMHB')]
    _bugs = _df[_df['Label_left'].str.contains('BUGS')]

    # bugmedia and bugbug control median values
    _bugs_media = _bugs[_bugs['Label_right'].str.contains('CAMHB')]
    bugs_media_median = _bugs_media[tp+'_norm2'].median()

    _bugs_bugs = _bugs[_bugs['Label_right'].str.contains('BUGS')]
    bugs_bugs_median = _bugs_bugs[tp+'_norm2'].median()

    # normalize these to bugs-media
    _abx_cp = _abx[_abx['Label_right'].str.contains('cp')]
    _abx_media = _abx[_abx['Label_right'].str.contains('CAMHB')]
    _bugs_cp = _bugs[_bugs['Label_right'].str.contains('cp')]
    bugmednorm = normalize(pd.concat([_abx_cp, _abx_media, _bugs_cp]),
    					   bugs_media_median)

    # normalize this to bugs-bugs
    _abx_abx = _abx[_abx['Label_right'].str.contains('abx')]
    _abx_bugs = _abx[_abx['Label_right'].str.contains('BUGS')]
    bugbugnorm = normalize(pd.concat([_abx_abx,_abx_bugs]), bugs_bugs_median)

    # do not normalize controls but add in relevant data
    _med_med = _med[_med['Label_right'].str.contains('CAMHB')]
    _med_med_median = _med_med[tp+'_norm2'].median()
    controlnorm = pd.concat([normalize(_bugs_bugs, bugs_bugs_median),
                            normalize(_bugs_media, bugs_media_median),
    						normalize(_med_med, _med_med_median)])

    _df = pd.concat([bugmednorm, bugbugnorm, controlnorm])

    return _df

def bootstrap_median(ar):
    ''' Resample from all dataframe rows '''
    return np.median(np.random.choice(ar, len(ar)))

def boot(df,boot_val,parallel=False):
    '''
    :df: df to bootstrap
    :boot_val: column name of values to bootstrap sample
    '''
    med_agg = bootstrap_median if parallel else 'median'
    grouped = df[['Label_left','Label_right',boot_val]].groupby(['Label_left','Label_right'])
    return grouped.aggregate(med_agg)

def boot_parallel(arglist):
    return boot(*arglist)

def SE_boot_(df):
    '''
    :df: normalized df
    :df: normalized df with se
    :MED: median normalized df with se
    '''
    df['comboID'] = df[['Label_left','Label_right']].apply(lambda x:\
    ','.join(x.dropna().astype(str)),axis=1)

    MED = df[['Label_left','Label_right','t0_norm2']].groupby(['Label_left',\
    'Label_right']).median()

    cores = mp.cpu_count()//int(POOL_NUM)
    pool = Pool(cores)
    boot_ = pool.map(boot_parallel,([df,'t0_norm2',j] for j in range(1,1001)))
    pool.close()
    pool.join()
    pool.clear()

    booted = pd.concat(boot_, axis=1)
    MED['t0_norm2_SE'] = np.std(booted.values,axis=1)
    MED = MED.reset_index()
    MED['comboID'] = MED[['Label_left','Label_right']].apply(lambda x:\
    ','.join(x.dropna().astype(str)), axis=1)
    SE_dict = dict(zip(MED.comboID, MED['t0_norm2_SE']))
    df['t0_norm2_SE'] = df['comboID'].map(SE_dict)
    return df,MED

##############################################################################
################################# NORMALIZE ##################################
qc_params = snakemake.input[2]
# get qc_params
qc_paramsdf = pd.read_csv(qc_params)
OUTLIER = qc_paramsdf.iloc[3,1]
POOL_NUM = qc_paramsdf.iloc[4,1]

if OUTLIER == 0:
    fil_file = snakemake.input[0]
elif OUTLIER == 1:
    fil_file = snakemake.input[1]

fil = pd.read_csv(fil_file,index_col=0)
day = fil_file.split('/')[-1].split('_')[0]
chip = fil_file.split('/')[-1].split('_qcfiltered')[0]
data_folder = 'data/output/'+day+'/'+chip
figure_folder='figures/'+day+'/'+chip

STRAIN = fil_file.split('/')[-1].split('_')[1]
CPP = fil_file.split('/')[-1].split('_')[2]

norm = get_growth_change(fil, tp='t0')

# bootstrap to get se of t0_norm2
norm,MED = SE_boot_(norm)

# add chip info
norm['chipID'] = chip

# save output
norm.to_csv(snakemake.output[0])
