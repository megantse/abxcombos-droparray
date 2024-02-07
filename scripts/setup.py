#!/usr/bin/env python3
'''
Script: setup.py
Last update: 2024 Feb, mtse

Description: SETUP CONFIG
    1. calc rescale vector for config
Input:
    1. data/raw/{DAY}/{exp}_cpp{n}/ # image directory
	2. data/raw/{DAY}/{DAY}_background/ # background image directory
	3. configs/{DAY}/template_config.yml # simplifies re-running from top for debugging
Output:
    1. figures/{DAY}/{exp}_cpp{n}_filter_test.png # sanity-check image
	2. configs/{DAY}/{exp}_cpp{n}_config.yml # new config
'''
##############################################################################
################################## IMPORTS ###################################
from skimage import io
import yaml
import sys, os, re, time
import argparse
import numpy as np
import matplotlib.pyplot as plt

##############################################################################
################################# FUNCTIONS ##################################
def get_max_drop(channel, q=0.99969):
	'''
	Find highest-intensity droplet in channel.

	Inputs:
	:channel: (NxN ndarray) pixel values from .tif file
	:q: (float) Quantile value between 0 and 1.
				0 = marks everything; 1 = marks nothing

	Outputs:
	:m: (float) Threshold pixel value
	'''

	ch = channel.flatten()
	m = np.quantile(ch, q)

	return m

def mark_max_drop(channel, q):
	'''
	Mark the selected quantile of pixels to check that
	they're in a droplet and not an artifact.
	Outputs intended for plotting.

	Inputs:
	:channel: (NxN ndarray) pixel values from .tif file
	:q: (float) Quantile value passed to get_max_drop.

	Outputs:
	:m: (float) get_max_drop output
	:cn: (NxN ndarray) binary-filtered channel; pix<m = 0; pix>m = 255
	:alphas: (NxN ndarray) Transparency values for plotting
	'''

	m = get_max_drop(channel, q)
	cn = np.uint64(channel > m) * 255
	alphas = cn/255

	return m, cn, alphas


def maxdrop_plots(rgb, q=[0.99969]*3, save=False):
	'''
	Plotting function for three-channel images.
	Default quantile value is approximately one droplet.

	Inputs:
	:rgb: (list) 3 NxN ndarrays corresponding to each channel
	:q: (list) List of quantile values to use for each channel
	:save: (str) directory to save image to. if False, does not save.

	Outputs: Max R G and B values.

	Plots images.
	'''

	if len(rgb)!=3:
		raise ValueError("Argument `rgb` must be a list of 3 numpy arrays")

	fig, axes = plt.subplots(1, 3, sharex='all', sharey='all', figsize=(30, 7))
	ax = axes.ravel()

	maxes = []

	for i in range(3):
		m, cn, alphas = mark_max_drop(rgb[i], q[i])
		maxes.append(m)
		ax[i].imshow(rgb[i], cmap='gray')
		ax[i].imshow(cn, cmap='Reds', alpha=alphas)
		ax[i].set_title('Channel '+str(i+1)+'\n'+str(100*q[i])+'th quantile')

	plt.show()

	if save:
		plt.savefig(save)

	return maxes

def rescale_cfg(cfg, maxes, rescale_info=['Reporter']):
	'''
	Updates config dictionary with rescale info and writes it to a new file.

	Inputs:
	:cfg: (dict) Old config dictionary
	:maxes: (list, length 3) Maximum value in each channel, from maxdrop_plots
	:cfg_dir: (str) Directory to write new config file to. Must end in .yml.
	:rescale_info: (list) Relevant additional info for rescale text.
						Loosely intended to provide info about any additional
						elements of the rescale vector.

	Modifies:
	:cfg: (dict) Same as before, but with 'rescale' and 'rescale text' added
				to it. Note that it is *not* a copy of the original.

	Writes new config file to :cfg_dir:
	'''

	mm = max(maxes)

	rescale = [float(round(mm/i, 1)) for i in maxes]
	rescale.append(1)

	cfg['image']['rescale'] = rescale
	cfg['image']['rescale_text'] = ['{:.2f}'.format(x) for x in maxes] + rescale_info

	return cfg

def path_cfg(cfg, img_dir, bg_dir):
	'''
	Creates new config file for one experiment from a template config.

	Inputs:
	:cfg: (dict) config dict, from snakemake
	:img_dir: (str) directory, relative to data/raw/, with experiment images
	:bg_dir: (str) directory, relative to data/raw/, with background images

	Outputs:
	:config: (dict) updated config dictionary

	'''
	config = cfg.copy()

	# below code for chip == 'data/raw/{DAY}/{exp}_cpp{n}/'
	# chip = img_dir[:-1].split('/')[-1]

	# below code for chip == first image file for this chip
	# 'data/raw/{{DAY}}/{{exp}}_cpp{{n}}/{{exp}}_cpp{{n}}{T}_{X}_{Y}.tif'
	chip = img_dir.split('/')[3]

	# change parameters:
	# path to image files
	config['image']['base_path'] = os.path.join(*(img_dir.split('/')[:4]))

	# path to background images
	config['image']['base_path2'] = bg_dir

	# get list of files in the image directory
	files = os.listdir(config['image']['base_path'])
	# print(files)

	# .tif names
	for tp in ['premerge', 't0']:
		# check that at least one image exists in the folder with that suffix
		config['image']['names'][tp] = chip + config['image']['names'][tp]

	return config

##############################################################################
################################### SETUP ####################################
if __name__ == '__main__':

	# make the new config file for the experiment
	cfg = path_cfg(snakemake.config, snakemake.input[2], snakemake.input[0])

	# read in file for rescaling calculation
	im = io.imread(os.path.join(
		cfg['image']['base_path'], 
		f"{cfg['image']['names']['premerge']}_4_3.tif"
	)) # why _4_3?
	rcol, gcol, bcol = im[...,0], im[...,1], im[...,2]

	mx = maxdrop_plots([rcol, gcol, bcol], save=snakemake.output[0])

	cfg = rescale_cfg(cfg, mx, ['Reporter', '6000'])

	cfg_dir='./configs/'+cfg['image']['base_path'][9:] + '_config.yml'

	# create new config file
	with open(cfg_dir, 'w') as file:
		yaml.dump(cfg, file)

	# give the file some time to save
	time.sleep(10)
