#!/usr/bin/env

# inputs
DAY = '20211229' # YYYYMMDD
BUGS = ['ab'] # all strains lower case, 1 species at a time
CPPNS = [1] # no double digit indexing
BUGTYPES = [1] # 1 or 2 for sensitve or resistant, respectively

################################################################################
# variables
KEYS_FILE = 'Keys-crl.xlsx'
NOTES_FILE = 'Notes-crl.xlsx'
SYNCTRL_FILE = 'synergy_control.xlsx'
QC_FILE = 'qc_params.csv'
QC_notebook = 0 # 0 for script; 1 for notebook
IMAGE_X_MAX = 5 # x-val of last .tif
IMAGE_Y_MAX = 8 # y-val of last .tif

# input dependent variables
CHIPS = expand("{day}_{bug}{bugtype}_cpp{n}", bug=BUGS, bugtype=BUGTYPES, n=CPPNS, day=DAY)
configfile: "configs/"+DAY+"/template_config.yml"
# BACKGROUND = list(range(2,config['image']['bg_num']+1))
################################################################################
rule all:
	input:
		expand('figures/{day}/{chip}_compiled_all_figures.pdf',chip=CHIPS,day=DAY) # 0
	params:
		day=DAY # 0
	shell:
		'bash scripts/vm_transfer_end.sh {params.day}'

# rule bg_avg:
# # change to output nd2 files in format compatible with rest of pipeline
# 	input:
# 		expand('data/raw/{{DAY}}/{{DAY}}_background/{bug}_'+config['image']['bg_input']+
# 				'{bg}.tif',bug=BUGS,bg=BACKGROUND), # 0 list of background images
# 	output:
# 		'data/raw/{DAY}/{DAY}_background/background_image.tif' # 0 averaged background image
# 	script:
# 		"scripts/background.py"

rule nd2_export:
	input:
		expand('data/raw/{{DAY}}/{{chip}}_{T}/', # 0: premerge image directory
			T=["premerge", "t0"]) 
	output:
		expand('data/raw/{{DAY}}/{{chip}}/{{chip}}{T}_{X}_{Y}.tif', 
			T=[config["image"]["names"]["premerge"], config["image"]["names"]["t0"]],
			X=range(1,IMAGE_X_MAX+1), Y=range(1,IMAGE_Y_MAX+1)),
		'data/raw/{DAY}/{chip}/nd2_metadata.pkl'
	script:
		'scripts/nd2_export.py'

rule nd2_background:
	input:
		'data/raw/{DAY}/{DAY}_background/'
	output:
		expand('data/raw/{{DAY}}/{{DAY}}_background/background_{X}_{Y}.tif', 
			X=range(1,IMAGE_X_MAX+1), Y=range(1,IMAGE_Y_MAX+1)),
		'data/raw/{DAY}/{DAY}_background/nd2_metadata.pkl'
	script:
		'scripts/nd2_export.py'

rule setup:
	input:
		'data/raw/{DAY}/{DAY}_background/', # 0 background image directory
		"configs/{DAY}/template_config.yml", # 1 simplifies re-running from top for debugging
		expand('data/raw/{{DAY}}/{{exp}}_cpp{{n}}/{{exp}}_cpp{{n}}{T}_{X}_{Y}.tif', 
			T=[config["image"]["names"]["premerge"], config["image"]["names"]["t0"]],
			X=range(1,IMAGE_X_MAX+1), Y=range(1,IMAGE_Y_MAX+1)),
		expand('data/raw/{{DAY}}/{{DAY}}_background/background_{X}_{Y}.tif', 
			X=range(1,IMAGE_X_MAX+1), Y=range(1,IMAGE_Y_MAX+1))
	output:
		'figures/{DAY}/{exp}_cpp{n}_filter_test.png', # 0 sanity-check image
		'configs/{DAY}/{exp}_cpp{n}_config.yml' # 1 new config
	script:
		"scripts/setup.py"

rule find_drops:
	input:
		'configs/{DAY}/{chip}_config.yml', # 0 new config file
		expand('data/raw/{{DAY}}/{{chip}}/{{chip}}{T}_{X}_{Y}.tif', 
			T=[config["image"]["names"]["premerge"], config["image"]["names"]["t0"]],
			X=range(1,IMAGE_X_MAX+1), Y=range(1,IMAGE_Y_MAX+1)),
		expand('data/raw/{{DAY}}/{{DAY}}_background/background_{X}_{Y}.tif', 
			X=range(1,IMAGE_X_MAX+1), Y=range(1,IMAGE_Y_MAX+1))
	output:
		'data/interim/{DAY}/{chip}_droplets_a.csv', # 0
		'data/interim/{DAY}/{chip}_droplets_b.csv', # 1
		'data/interim/{DAY}/{chip}_droplets_c.csv', # 2
		'data/interim/{DAY}/{chip}_droplets_d.csv' # 3 droplets csv
	script:
		"scripts/find_drops.py"

rule make_cluster:
	input:
		'data/interim/{DAY}/{chip}_droplets_d.csv', # 0
		'notebooks/master_interactive_cluster.ipynb', # 1
		expand('data/raw/{{DAY}}/{{chip}}/{{chip}}{T}_{X}_{Y}.tif', 
			T=[config["image"]["names"]["premerge"], config["image"]["names"]["t0"]],
			X=range(1,IMAGE_X_MAX+1), Y=range(1,IMAGE_Y_MAX+1))
	output:
		'notebooks/{DAY}/{chip}_interactive_cluster.ipynb'
	shell:
		'cp notebooks/master_interactive_cluster.ipynb {output}'

rule cluster:
	input:
		expand('data/interim/{{DAY}}/{chip}_droplets_d.csv',chip=CHIPS), # 0
		expand('notebooks/{{DAY}}/{chip}_interactive_cluster.ipynb',chip=CHIPS), # 1
		expand('configs/{{DAY}}/{chip}_config.yml',chip=CHIPS), # 2
		NOTES_FILE # 3
	params:
		day=DAY, # 0
		chips=CHIPS # 1
	output:
		'data/interim/{DAY}/cluster_notebook_output.txt', # 0
		expand('data/interim/{{DAY}}/{chip}_clustered.csv',chip=CHIPS), # 1
		expand('data/interim/{{DAY}}/{chip}_missing_clusters.csv',chip=CHIPS) # 2
	shell:
		'jupyter notebook --ip=0.0.0.0 --port=8888 --notebook-dir="notebooks" '
		'--no-browser --allow-root; '
		'bash scripts/check_cluster.sh {params.day} {params.chips}'

rule pre_post:
	input:
		'data/interim/{DAY}/{chip}_clustered.csv', # 0
		'configs/{DAY}/{chip}_config.yml', # 1
		'data/interim/{DAY}/cluster_notebook_output.txt' # 2
	output:
		'data/interim/{DAY}/{chip}_pre_post.csv', # 0
		'data/interim/{DAY}/{chip}_condensed.csv' # 1
	script:
		"scripts/pre_post.py"

rule qc_preprocess:
	input:
		'data/interim/{DAY}/{chip}_clustered.csv', # 0
		'data/interim/{DAY}/{chip}_condensed.csv' # 1
	output:
		'data/interim/{DAY}/{chip}_trimmed.csv' # 0
	script:
		"scripts/qc_preprocess.py"

rule make_qc:
	input:
		'data/interim/{DAY}/{chip}_trimmed.csv', # 0
		'notebooks/master_qc_plots.ipynb' # 1
	output:
		'notebooks/{DAY}/{chip}_qc_plots.ipynb' # 0
	shell:
		'cp notebooks/master_qc_plots.ipynb {output}'

if QC_notebook == 0:
	rule qc_filter:
		input:
			'data/interim/{DAY}/{chip}_trimmed.csv', # 0
			'data/interim/{DAY}/{chip}_pre_post.csv', # 1
			QC_FILE # 2
		output:
			'data/output/{DAY}/{chip}_qcfiltered.csv' # 0
		script:
			"scripts/qc_filter.py"

elif QC_notebook == 1:
	rule qc_filter:
		input:
			expand('data/interim/{{DAY}}/{chip}_trimmed.csv',chip=CHIPS), # 0
			expand('data/interim/{{DAY}}/{chip}_pre_post.csv',chip=CHIPS), # 1
			expand('notebooks/{{DAY}}/{chip}_qc_plots.ipynb',chip=CHIPS) # 2
		params:
			day=DAY, # 0
			chips=CHIPS # 1
		output:
			"data/interim/{DAY}/qc_notebook_output.txt", # 0 dummy file
			expand('data/output/{{DAY}}/{chip}_qcfiltered.csv',chip=CHIPS) # 1
		shell:
			'jupyter notebook --ip=0.0.0.0 --port=8888 --notebook-dir="notebooks" '
			'--no-browser --allow-root; bash scripts/check_qc.sh {params.day} {params.chips}'

rule outlier_filter:
	input:
		'data/output/{DAY}/{chip}_qcfiltered.csv' # 0
	output:
		'data/output/{DAY}/{chip}_outlierfiltered.csv' # 0
	script:
		"scripts/outlier_filter.py"

rule normalize:
	input:
		'data/output/{DAY}/{chip}_qcfiltered.csv', # 0
		'data/output/{DAY}/{chip}_outlierfiltered.csv', # 1
		QC_FILE # 2
	output:
		'data/output/{DAY}/{chip}_normalized.csv' # 0
	script:
		"scripts/normalize.py"

rule abx_only:
	input:
		'data/output/{DAY}/{chip}_normalized.csv', # 0
		KEYS_FILE, # 1
		QC_FILE # 2
	output:
		'data/output/{DAY}/{chip}_abxonly_droplet_data.csv', # 0
		'data/output/{DAY}/{chip}_abxonly_summary_data.csv' # 1
	script:
		"scripts/abx_only.py"

rule cp_only:
	input:
		'data/output/{DAY}/{chip}_normalized.csv', # 0
		KEYS_FILE, # 1
		QC_FILE # 2
	output:
		'data/output/{DAY}/{chip}_cponly_droplet_data.csv', # 0
		'data/output/{DAY}/{chip}_cponly_summary_data.csv', # 1
		'data/output/{DAY}/{chip}_cponly_high_activity.csv' # 2
	script:
		"scripts/cp_only.py"

rule abx_abx:
	input:
		'data/output/{DAY}/{chip}_normalized.csv', # 0
		KEYS_FILE, # 1
		SYNCTRL_FILE, # 2
		QC_FILE # 3
	output:
		'data/output/{DAY}/{chip}_abxabx_droplet_data.csv', # 0
		'data/output/{DAY}/{chip}_abxabx_summary_data.csv', # 1
		'data/output/{DAY}/{chip}_abxabx_negctrl_droplet_data.csv', # 2
		'data/output/{DAY}/{chip}_abxabx_negctrl_summary_data.csv', # 3
		'data/output/{DAY}/{chip}_abxabx_negctrl_summed_summary_data.csv', # 4
		'data/output/{DAY}/{chip}_abxabx_summed_summary_data.csv', # 5
		'data/output/{DAY}/{chip}_abxabx_synctrl.csv' # 6
	script:
		"scripts/abx_abx.py"

rule abx_cp:
	input:
		'data/output/{DAY}/{chip}_normalized.csv', # 0
		'data/output/{DAY}/{chip}_abxabx_synctrl.csv', # 1
		KEYS_FILE, # 2
		QC_FILE # 3
	output:
		'data/output/{DAY}/{chip}_abxcp_droplet_data.csv', # 0
		'data/output/{DAY}/{chip}_abxcp_negctrl_droplet_data.csv', # 1
		'data/output/{DAY}/{chip}_abxcp_negctrl_summary_data.csv', # 2
		'data/output/{DAY}/{chip}_abxcp_negctrl_summed_summary_data.csv', # 3
		'data/output/{DAY}/{chip}_abxcp_summed_summary_data.csv', # 4
		'data/output/{DAY}/{chip}_abxcp_summary_data.csv', # 5
		'data/output/{DAY}/{chip}_abxcp_chip_hit.csv', # 6
		'figures/{DAY}/{chip}_abxcp_chiphits.pdf' # 7
	script:
		"scripts/abx_cp.py"

rule stats:
	input:
		'data/output/{DAY}/{chip}_normalized.csv', # 0
		 KEYS_FILE, # 1
		'data/interim/{DAY}/{chip}_missing_clusters.csv', # 2
		'data/interim/{DAY}/{chip}_droplets_b.csv', # 3
		'data/interim/{DAY}/{chip}_clustered.csv', # 4
		'data/interim/{DAY}/{chip}_trimmed.csv', # 5
		'data/output/{DAY}/{chip}_qcfiltered.csv', # 6
		'data/output/{DAY}/{chip}_outlierfiltered.csv', # 7
		'data/interim/{DAY}/{chip}_pre_post.csv', # 8
		'data/output/{DAY}/{chip}_cponly_summary_data.csv', # 9
		'data/output/{DAY}/{chip}_abxonly_summary_data.csv', # 10
		'data/output/{DAY}/{chip}_abxabx_summed_summary_data.csv', # 11
		'figures/{DAY}/{chip}_abxcp_chiphits.pdf' # 12
	output:
		'data/output/{DAY}/{chip}_stat_zprime_st.csv', # 0
		'data/output/{DAY}/{chip}_stat_zprime_abx.csv', # 1
		'figures/{DAY}/{chip}_compiled_all_figures.pdf' # 2
	script:
		"scripts/stats.py"
