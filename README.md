How to analyze DropArray screening data

Packages:
- listed in snakemake.yml file


Example below is using files copied from a storage unit called bucket and running this analysis using a virtual machine.


Step 0: Prerequisites (skip, when copying ready VM)

In order to run the pipeline, you need a dedicated folder in your VM. Mine is set to ~/analysis/, though you can call it whatever works for you. I will be referring to that folder throughout; all commands below are listed as though you’re running them from that directory. You only need to do this part once.


Organize this directory as follows:

analysis/

├── ./Keys-CRL.xlsx

├── ./Notes-CRL.xlsx

├── ./crl_props.xlsx

├── ./Snakefile

├── ./configs

│   └── ./configs/template_config.yml

├── ./dag.pdf

├── ./data/

│   ├── ./data/interim

│   ├── ./data/output

│   └── ./data/raw

├── ./figures/

├── ./kchip_py3

│   └── [scripts for kchip]

├── ./logs/

├── ./notebooks/

│   ├── ./notebooks/master_interactive_cluster.ipynb

│   └── ./notebooks/master_qc_plots.ipynb

├── ./qc_params.csv

└── ./scripts/

    └── [scripts for antibiotic-compound combos analysis]
    

Create conda environment

$ conda env create -f snakemake.yml


Step 1: Preparation

Copy files to directory. 


First make a new target directory.

$ mkdir data/raw/[day]/


Mount your bucket.

$ gcsfuse --implicit-dirs production-u19 ~/bucket


Generally I like to do each day in chunks by bug species since the VM has limited storage capacity. You can usually store about two bugs’ worth of images (for the sensitive and resistant strains of a single species) in the folder before you start having issues. These commands usually take 10-15 minutes to run depending on how much you’re moving.

$ cp -r ~/bucket/data/raw/[day]/[day]_[bug]_[cpp#] data/raw/[day]/


Also be sure to copy the background images for the day if you haven’t already:

$ cp -r ~/bucket/data/raw/[day]/[day]_background data/raw/[day]/

Set up config file

Copy template config to folder for the new day:

$ mkdir configs/[day]

$ cp configs/template_config.yml configs/[day]/



Open the new copy of the file and adjust the values as necessary using Vim.  

$ vi configs/[day]/template_config.yml

$ i 

Modify the day number, bacteria abbreviation, and strain number

*Press control+c or Exc to get out of insert mode*

$ :wq # to save changes

$ :q! # to not save changes


Usually, this will mean changing the values in image/names/. The default values in the template are the most commonly found name infixes. In other words, if the listed premerge value is _premerge, then we anticipate premerge files with name pattern: 

[day]_[bug]{1 or 2}_cpp[num]_premerge_[X]_[Y].tif

[day]_[bug]{1 or 2}_cpp[num]_tfinal_[X]_[Y].tif

The pipeline will autofill the day, bug, chip, and X/Y identifiers. If the infix (bolded portion above) is different from the value in the template file, change it to match. Note that only one underscore is bolded here. 


Similarly, adjust the bg_input value. eg. bg_input: background implies background{1-10}.tif input images. 1to1_bkg: background implies one-to-one background subtraction using corresponding background image tile. 

Set up Snakefile

Open the Snakefile. Edit the first three lines in order to get it to take in the correct data: 

DAY = '20211229' 

BUGS = ['ab']

CPPNS = [1,2]

BUGTYPES = [1, 2]

The values above tell Snakemake to run on the following raw datasets:

20211229_ab1_cpp1, 20211229_ab1_cpp2, 20211229_ab2_cpp1, 20211229_ab2_cpp2


Note that BUGS is a list with only one value. You can run more than one bug by adding a second bug identifier, but I wouldn’t run on more than 8 datasets at once, depending on how much space you have available on your VM. If you run out of space mid-analysis, the pipeline will auto-exit, which can be a real pain to debug. 

Also note that the pipeline will construct all possible name combinations as shown above. If any of them doesn’t exist, it will error out. 


Step 2: Running the pipeline

Continue to do everything in the ~/analysis/pilot_screen directory FYI, is this correct?



Activate the environment

$ conda activate snakemake


Final check (dry run and port check)

Snakemake has a function that will check to make sure your pipeline will run the way you want it to, without touching any of the files yet. 

$ snakemake -n

You should see a list of jobs at the end of the command output.


Start the pipeline

$ snakemake --cores 32

Use --unlock if needed

The pipeline will immediately start. 


Clustering: When it reaches the clustering step, it will create a clustering notebook for each dataset and export them to port 8888 on your machine. Open the clustering notebooks like you would any other, navigate to the appropriate day, and run through clustering for each one. Once satisfied with the results, quit the Jupyter session and the pipeline will continue running. If it doesn’t detect the right output files, the pipeline will stop and shut down. 


QC: 

IF NOTEBOOK==1:

At the end, it will create a quality-control notebook for each dataset. This follows the same process as the clustering step above--if the appropriate output files aren’t detected by the end, it will produce an error.


IF NOTEBOOK==0:

Run through qc_filter.py scripts based on parameters in qc_params.csv.

Step 3: Closing the pipeline

Transferring output data files and figures to bucket

Rule all w/ script “vm_transfer_end.sh” will copy everything to the bucket and shutdown the VM.

To do it individually or without vm_transfer_end.sh:


$ cp -r ~/analysis/data/output/[day]/*.csv ~/bucket/data/output

$ cp -r ~/analysis/data/interim/[day]/*.csv ~/bucket/data/interim

$ cp -r ~/analysis/configs/[day]/*.yml ~/bucket/configs

$ cp -r ~/analysis/notebooks/[day]/*.ipynb ~/bucket/notebooks

$ cp -r ~/analysis/figures/[day]/*.png ~/bucket/figures



Make sure to delete the input files from the VM so you can load in new ones and re-run.


$ rm -r ~/analysis/data/*/[day]/

$ rm -r ~/analysis/*/[day]/
