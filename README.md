---
Title: "Replaying the evolutionary tape to investigate subgenome dominance in allopolyploid _Brassica napus_"
Methylome Analysis: "Chad Niederhuth"
Collaborators: "Kevin A. Bird, Shujun Ou, Malia Gehan, J. Chris Pires, Zhiyong Xiong, Robert VanBuren Patrick P. Edger"
Raw data: "Link to SRA data will be posted at later date"
---
This repository is for scripts and processed data for the paper:

Please cite this paper if you use any of the resources here.  

All analyses performed on the Michigan State University High Performance Computing Cluster (HPCC)

To reproduce the analysis, follow these steps:

**NOTE #1:** This analysis assumes you will be using Anaconda and I have provided a yml file to easily create an environment for repeating analysies. 

1) Clone this git repository

git clone https://github.com/niederhuth/The-molecular-basis-of-kale-domestication-Comparative-transcriptomics
cd Bnapus-polyploidy

2) Create the conda environment

conda env create -f Bnapus-polyploidy.yml

3) You will now need to create a symbolic link within this environment for methylpy to work. 

cd /env/Bnapus-polyploidy/lib
ln -s libgsl.so.23.0.0 libgsl.so.0

4) Return to the cloned git repository

cd ~/Bnapus-polyploidy

5) Create a data folder and cd to it

mkdir data
cd data

6) Run the setup.sh script (see note #2)

bash ../scripts/setup.sh

7) For each sample, run methylpy

cd <sample_directory>
bash ../../scripts/run_methylpy.sh

or submit as a job

8) When methylpy is finished, for each sample, you can run the various scripts in the "analyses_sh" directory (see note #2)

cd <sample_directory>
bash ../../scripts/analyses_sh/<script_name>

or submit as a job

These will create a series of outputs for each analysis in <sample_name>/combined/results

These should match data found in the figures_tables directory

9) You can run the methylation_analysis.R in the figures_tables directory to output plots and other analyzed results

cd figures_tables
Rscript ../scripts/methylation_analysis.R

**NOTE #2:** These scripts were written for use on the MSU HPCC. To run them on your computer or a different environment, you will need to change the header of each shell script (python and R scripts do not require changing) to something appropriate for your system. In each shell script, you will also need to either modify or delete these lines to a path appropriate for your system:
export PATH="$HOME/miniconda3/envs/Bnapus-polyploidy/bin:$PATH"
export LD_LIBRARY_PATH=""$HOME/miniconda3/envs/Bnapus-polyploidy/lib:$LD_LIBRARY_PATH"

You also may need to delete this line

cd $PBS_O_WORKDIR
