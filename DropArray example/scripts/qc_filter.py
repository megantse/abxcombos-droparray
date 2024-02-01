'''
Script: qc_filter.py
Last update: 2021 March, mzhu

Description: QC FILTERING
    1. filter on droplet area size
    2. filter on barcode distance
Input:
    1. data/interim/{DAY}/{chip}_trimmed.csv
    2. data/interim/{DAY}/{chip}_pre_post.csv
    3. qc_params.csv
Output:
    1. data/output/{DAY}/{chip}_qcfiltered.csv
'''
##############################################################################
################################## IMPORTS ###################################
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
from plot_param import parameters
parameters()

##############################################################################
################################ QC FILTER ###################################
trimmed_file = snakemake.input[0]
prepost_file = snakemake.input[1]
qc_params = snakemake.input[2]

# get qc_params
qc_paramsdf = pd.read_csv(qc_params)
Area_lowerbound = qc_paramsdf.iloc[0,1]
Area_upperbound = qc_paramsdf.iloc[1,1]
Barcode_dist = qc_paramsdf.iloc[2,1]

trimmed = pd.read_csv(trimmed_file,index_col=0)
pre_post = pd.read_csv(prepost_file,index_col=0)
day = trimmed_file.split('/')[-1].split('_')[0]
chip = trimmed_file.split('/')[-1].split('_trimmed')[0]
data_folder = 'data/output/'+day+'/'+chip
figure_folder='figures/'+day+'/'+chip

# plot histogram of droplet area size
plt.figure()
plt.hist(trimmed['t0_Area'],bins=50,range=(0,1000), label = 't0', alpha=0.4)
plt.xlabel('Area')
plt.ylabel('Frequency')
plt.title('Droplet area: untrimmed')
plt.xlim(200,800)
plt.legend()
plt.savefig(figure_folder+'_droplet_area.png',bbox_inches='tight',dpi=300)

# plot barcode distances
Ldist = trimmed[['Label_left','distance_left']]
Rdist = trimmed[['Label_right','distance_right']]

f, axes = plt.subplots(1, 2, figsize=(40,5))
sns.boxplot(x='Label_left',y='distance_left',data=Ldist, ax=axes[1])
sns.boxplot(x='Label_right',y='distance_right',data=Rdist, ax=axes[0])
for ax in f.axes:
    ax.tick_params(labelrotation=90)
ylim = [0.0,0.05]
plt.setp(axes, ylim=ylim)
f.suptitle('Barcode distances: untrimmed')
plt.savefig(figure_folder+'_barcode_distances.png',bbox_inches='tight',dpi=300)

# filter based on barcode distance
Ldistance_trimmed = trimmed[trimmed['distance_left']<Barcode_dist]
distance_trimmed = Ldistance_trimmed[Ldistance_trimmed['distance_right']<Barcode_dist]

# plot filtered clusters
plt.figure()
plt.scatter(distance_trimmed['PlaneX_left'],distance_trimmed['PlaneY_left'],alpha=0.05, label = 'label_left')
plt.scatter(distance_trimmed['PlaneX_right'],distance_trimmed['PlaneY_right'],alpha=0.05, label = 'label_right')
plt.gca().set_title('On_plane coordinates of filtered droplets from wells')
plt.legend()
plt.savefig(figure_folder+'_barcode_distance_trimmed_clusters.png',bbox_inches='tight',dpi=300)

# filter based on droplet area size
d_area_trimmed = distance_trimmed[distance_trimmed['t0_Area']<Area_upperbound]
d_area_trimmed = d_area_trimmed[d_area_trimmed['t0_Area']>Area_lowerbound]
d_area_trimmed['chipID'] = chip

# save output
d_area_trimmed.to_csv(snakemake.output[0])

# plot locations of all identified droplets on 2D plane to represent the chip
# can identify portion of chip where droplet identification failed
# regions may due to cross merging
fig, axes = plt.subplots(2,1,figsize=(10,20))

axes[0].plot(pre_post['Pre_GlobalX'],-pre_post['Pre_GlobalY'],'.',ms=1)
axes[0].set_title('Wells Detected')

axes[1].plot(pre_post['Post_GlobalX'],-pre_post['Post_GlobalY'],'.',ms=1)
axes[1].set_title('Wells Detected - Edge Wells Removed')

plt.tight_layout()
plt.savefig(figure_folder+'_droplet_locations.png',bbox_inches='tight',dpi=300)

# plot label counts post filtering
left_count = d_area_trimmed['Label_left'].value_counts()
right_count = d_area_trimmed['Label_right'].value_counts()
counts = left_count + right_count

fig, axes = plt.subplots(1,1,figsize=(50,4))
axes.plot(counts.values,'o')
axes.set_xticks(range(len(counts.index.values)))
axes.set_xticklabels(counts.index.values,size=10,rotation=90)
plt.ylim(0, 1.5*max(counts))
axes.set_title('Counts - postfilter',size=20)
plt.savefig(figure_folder+'_cluster_counts_postfilter.png',bbox_inches='tight',dpi=300)

# plot control signal distributions
bug_drops = d_area_trimmed[(d_area_trimmed.Label_left.str.contains('BUGS')) & \
(d_area_trimmed.Label_right.str.contains('BUGS'))]
bugmed_drops = d_area_trimmed[(d_area_trimmed.Label_left.str.contains('BUGS')) & \
(d_area_trimmed.Label_right.str.contains('CAMHB'))]
media_drops = d_area_trimmed[(d_area_trimmed.Label_left.str.contains('CAMHB')) & \
(d_area_trimmed.Label_right.str.contains('CAMHB'))]
plt.figure()
plt.hist(bug_drops.t0_norm2,label='BUGS')
plt.hist(bugmed_drops.t0_norm2,label='BUG+CAMHB')
plt.hist(media_drops.t0_norm2,label='CAMHB')
plt.title('Controls')
plt.legend(bbox_to_anchor=(1.25,1))
plt.savefig(figure_folder+'_control_distributions.png',bbox_inches='tight',dpi=300)
