'''
Script: background.py
Last update: 2021 March, mzhu

Description: BACKGROUND IMAGE AVERAGING
    - Averages across multiple images (works for single images or stacks)
    - Eg. to average foreground images
    - Select settings
    - Sigma is for the gaussian filter
    - Use tuple for image stacks:
        - A 2048x2048 RYGBU imports as (2048, 2048, 5)
    - Use eg. sigma = (50,50,0) to filter in X and Y in each image, but not across channels
    - For single images, use eg. sigma = (50,50) to filter in X and Y
Input:
    1. expand('data/raw/{{DAY}}/{{DAY}}_background/'+config['image']['bg_input']+'{bg}.tif',bg=BACKGROUND)
Output:
    1. path.join(config['image']['base_path2'],config['image']['background_image']+'.tif')
'''
##############################################################################
################################## IMPORTS ###################################
import numpy as np
from skimage import io
from os import path
import scipy.ndimage as ndimage
import tifffile
import yaml
import warnings
warnings.filterwarnings('ignore')

##############################################################################
####################### BACKGROUND AVERAGING #################################
# sigma for gaussian filter
sigma = (10,10,0)

# create list to store image arrays in
images = []

# read in images
bkg_num_images = 10
input_img_names = snakemake.input
print(input_img_names)
# import each image and apply gaussian filter
# add the smoothed image to the list
for img in input_img_names:
    img_ = io.imread(img)#[:4,:,:]
    img_ = ndimage.gaussian_filter(img_, sigma=sigma, order=0)
    images.append(img_)

# stack images from list to ndarray
img_stack = np.stack(images)

# average the image values across the images
img_ave = img_stack.mean(axis=0)

# convert the ndarray data to a convenient output data type and shape
img16 = img_ave.astype('uint16')
if img16.shape[1] > img16.shape[2]:
    img16trans = np.transpose(img16, (2,0,1))
else:
    img16trans = img16

new_name = snakemake.output[0]

if len(img16.shape) == 3:
    tifffile.imsave(new_name, img16trans, photometric='minisblack')
else:
    tifffile.imsave(new_name, img16)
