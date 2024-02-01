from pims import ND2_Reader
from nd2reader import ND2Reader
import skimage.io as io
import os
import shutil
import pandas as pd
import math as mt
import numpy as np
import argparse
import time
import glob
import warnings
warnings.filterwarnings('ignore')



ROWS = 7
COLS = 10

# Exports nd2 file into individual tif filters
# python export_nd2_multi.py --day DAY --bug BUG --cppn cpp1_cpp2_cpp3

# def get_args():
#     parser = argparse.ArgumentParser()
#     parser.add_argument("--day",
#                         help="date of screen",
#                         type=str)
#     parser.add_argument("--bug",
#                         help="ab, kp, or pa",
#                         type=str)
#     parser.add_argument("--cppn",
#                         help="e.g. cpp1_cpp2_cpp3",
#                         type=str)
#     args = parser.parse_args()
#     return args

# import multiprocessing as mp
# from pathos.multiprocessing import ProcessingPool as Pool

def get_meta(fov):
    return list(fov.metadata.values())

def metadata_df(nd2):
    # create the multiprocessing pool
    cores = mp.cpu_count()//int(1)
    pool = Pool(cores)
    # process the dataFrame by mapping function to each df across the pool
    meta_list = pool.map(get_meta,[nd2[i] for i in range(nd2.sizes['m'])])
    # close down the pool and join
    pool.close()
    pool.join()
    pool.clear()

    meta = pd.DataFrame(meta_list,columns=nd2[0].metadata.keys())
    meta['delta_y'] = np.concatenate((np.asarray([0]),np.diff(meta['y_um'])))
    # rows = mt.ceil(meta['delta_y'].sum()/meta['delta_y'].max())
    meta['deltanorm_y'] = abs(meta['delta_y']/meta['y_um']*100)
    rows = meta[meta['deltanorm_y']>10].shape[0]+1
    cols = int(nd2.sizes['m']/rows)
    print(rows,cols)
    return rows, cols, meta

def save_img(vars_):
    m_, file_name = vars_
    nd2 = ND2_Reader(nd2file_global)
    nd2.iter_axes = 'm'
    nd2.bundle_axes = ['c','y','x']
    file = os.path.join(savepath_global,file_name)
    io.imsave(file,nd2[m_])
    return file

def export_background(day):
    """Export background images to TIF if they haven't been already."""

    savepath = os.path.join('./data/raw/',day,day+'_background')

    background_tifs = glob.glob(
        os.path.join(
            savepath, '*.tif'
        ))

    background_nd2 = glob.glob(os.path.join(savepath,'*','*'+'.nd2'))
    process_nd2_py(background_nd2[0], 'background', savepath)



def get_nd2_filepaths(day, bug, cpp, timepoint):
    base_name = f'{day}_background' if timepoint == 'background' else f'{day}_{bug}_{cpp}_{timepoint}'
    return glob.glob(
        os.path.join(
            './data/raw/',day,
            base_name,
            '*','*'+'.nd2'
            )
        )[0]

def process_nd2_pims(nd2file, root_name):
    """***DOES NOT WORK*** Export ND2 files as TIFs.

    One ND2 --> 40 TIFs."""
    raise NotImplementedError

    """with ND2_Reader(nd2file) as nd2:
                    print('created ND2 reader object')
                    nd2.iter_axes = 'm'
                    nd2.bundle_axes = ['c','y','x']

                    meta = pd.DataFrame(index=np.arange(0,len(nd2)))
                    print('created base dataframe')
                    meta['root_name'] = root_name+'_'
                    meta['m'] = meta.index
                    meta['tile_id'] = meta.m+1
                    meta['tile_row'] = meta.tile_id.apply(lambda x: int(mt.ceil(float(x)/COLS)))
                    meta['tile_col'] = meta.tile_id.apply(lambda x: int(mt.ceil(float(x)-((mt.ceil(float(x) / COLS)-1)*COLS))))
                    meta['tile_row_odd'] = meta['tile_row'] % 2
                    meta['tile_col_new'] = np.zeros
                    meta['tile_col_new'][meta['tile_row_odd']==1] = meta['tile_col']
                    meta['tile_col_new'][meta['tile_row_odd']==0] = COLS+1-meta['tile_col']
                    meta['rect_filename'] = (
                        meta['root_name']+meta['tile_row'].astype(str)+'_'+meta['tile_col_new'].astype(str)+'.tif'
                    )
                    print(len(meta['rect_filename']))
                    meta['enum_filename'] = meta['root_name']+'#'+meta['tile_id'].astype(str)+'.tif'
                    meta.to_pickle(os.path.join(savepath,'nd2_metadata.pkl'))

                    print('saving ',savepath_global,' ',root_name)
                    for j, fov in enumerate(nd2[:]): # this loads the entire file into memory
                        # remove [:] in line above and it should iterate over m-axis
                        filename = meta.loc[j+1,'rect_filename']
                        io.imsave(os.path.join(savepath_global,filename),fov)"""

def process_nd2_py(nd2file, root_name, savepath):
    """Export ND2 files as TIFs.

    One ND2 --> 40 TIFs."""
    print(nd2file)
    with ND2Reader(nd2file) as nd2:
        # switch 'm' with 'v'
        nd2.iter_axes = 'v'
        nd2.bundle_axes = 'cyx'

        meta = pd.DataFrame(index=np.arange(0,len(nd2)))
        print('created base dataframe')
        meta['root_name'] = root_name+'_'
        meta['m'] = meta.index
        meta['tile_id'] = meta.m+1
        meta['tile_row'] = meta.tile_id.apply(lambda x: int(mt.ceil(float(x)/COLS)))
        meta['tile_col'] = meta.tile_id.apply(lambda x: int(mt.ceil(float(x)-((mt.ceil(float(x) / COLS)-1)*COLS))))
        meta['tile_row_odd'] = meta['tile_row'] % 2
        meta['tile_col_new'] = np.zeros
        meta['tile_col_new'][meta['tile_row_odd']==1] = meta['tile_col']
        meta['tile_col_new'][meta['tile_row_odd']==0] = COLS+1-meta['tile_col']
        meta['rect_filename'] = (
            meta['root_name']+meta['tile_row'].astype(str)+'_'+meta['tile_col_new'].astype(str)+'.tif'
        )
        print(len(meta['rect_filename']))
        meta['enum_filename'] = meta['root_name']+'#'+meta['tile_id'].astype(str)+'.tif'
        meta.to_pickle(os.path.join(savepath,'nd2_metadata.pkl'))

        print('saving ',savepath,' ',root_name)
        for j, fov in enumerate(nd2[:]): # this loads the entire file into memory
            # remove [:] in line above and it should iterate over m-axis
            filename = meta.loc[j,'rect_filename']
            io.imsave(os.path.join(savepath,filename),fov)

def wait_for_files(input_args):
    """Wait for files to complete saving to appease Snakemake.

    Args:
    :input_args: list of arguments.
        if background: [date, 'background']
        else: [date, bug, cpp, timepoint]
    """

    if len(input_args) == 2: # background
        savepath = "_".join(input_args)
        target_num_files = ROWS*COLS
    else: # image
        savepath = "_".join(input_args[:3])
        target_num_files = ROWS*COLS*2 # two timepoints

    searchpath = os.path.join('./data/raw', input_args[0], savepath, '*.tif')

    i = 0
    while i < target_num_files:
        time.sleep(10)
        i = len(glob.glob(searchpath))

def delete_nd2():
    """Delete nd2 file directories after exporting to TIF."""
    for dir in snakemake.input:
        shutil.rmtree(os.path.join('.', dir))

def main(day,bug,cpp,cfg):

    timepoints = [cfg['image']['names']['premerge'], cfg['image']['names']['t0']]
    nd2_files = [get_nd2_filepaths(day, bug, cpp, tp) for tp in ['premerge','t0']]
    root_names = [f'{day}_{bug}_{cpp}{tp}' for tp in timepoints]

    print(nd2_files)
    print(root_names)

    for nd2file, root_name in zip(nd2_files, root_names):
        savepath = os.path.join('./data/raw/',day,day+'_'+bug+'_'+cpp)

        os.makedirs(savepath,exist_ok=True)

        process_nd2_py(nd2file, root_name, savepath)

            # # create the multiprocessing pool
            # cores = mp.cpu_count()//int(1)
            # pool = Pool(cores)
            # # process the dataFrame by mapping function to each df across the pool
            # img_list = pool.map(save_img,list(zip(meta['m'].values,meta['rect_filename'].values)))
            # # close down the pool and join
            # pool.close()
            # pool.join()
            # pool.clear()

if __name__=="__main__":

    start_time = time.time()
    # args = get_args().__dict__
    # main(**args)



    input_args = snakemake.input[0].split('/')[-2].split('_')
    if 'background' in input_args:
        export_background(input_args[0])
    else:
        day, bug, cpp, _ = input_args
        main(day, bug, cpp, snakemake.config)

    print("--- %s seconds ---" % (time.time() - start_time))

    print("Waiting for files to save...")

    wait_for_files(input_args)

    print('Done saving')
    if 'background' in input_args:
        pass
    else:
        delete_nd2()
