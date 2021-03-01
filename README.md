---
Title: "Replaying the evolutionary tape to investigate subgenome dominance in allopolyploid Brassica napus"
Authors: "Kevin A. Bird, Chad E. Niederhuth, Shujun Ou, Malia Gehan, J. Chris Pires, Zhiyong Xiong, Robert VanBuren Patrick P. Edger"
Raw data: "SRA Accession: PRJNA577908"
---
This repository is for DNA methylation scripts and data for the paper: https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.17137

Please cite this paper if you use any of the resources here. 

All analyses performed on the Michigan State University High Performance Computing Cluster (HPCC)

To reproduce the analysis, follow these steps:

**NOTE #1:** This analysis assumes you will be using Anaconda and I have provided a yml file to easily create an environment for repeating analyses. 

**1)** Clone this git repository

```
git clone https://github.com/niederhuth/Replaying-the-evolutionary-tape-to-investigate-subgenome-dominance  
cd Replaying-the-evolutionary-tape-to-investigate-subgenome-dominance   
```

**2)** Create the conda environment

```
conda env create -f scripts/Bnapus-polyploidy.yml
```

**3)** You will now need to create a symbolic link within this environment for methylpy to work. This will require you to cd into the environment located
in your anaconda (or miniconda) directory. 

```
cd miniconda3/env/Bnapus-polyploidy/lib  
ln -s libgsl.so.23.0.0 libgsl.so.0  
```

**4)** Return to the cloned git repository

```
cd ~/Replaying-the-evolutionary-tape-to-investigate-subgenome-dominance  
```

**5)** Create a data folder and cd to it

```
mkdir data
cd data
```

**6)** Run the setup.sh script (see note #2)

```
bash ../scripts/sh/setup.sh
```

**7)** For each sample, run methylpy

```
cd <sample_directory>  
bash ../../scripts/sh/run_methylpy.sh
```

or submit as a job

**8)** When methylpy is finished, for each sample, you can run the various scripts in the "sh" directory (see note #2)

```
cd **<sample_directory>**  
bash ../../scripts/_sh/**<script_name>**  
```

or submit as a job

These will create a series of outputs for each analysis in __**<sample_name>**/combined/results__

These should match data found in the figures_tables directory

**9)** You can run the methylation_analysis.R in the figures_tables directory to output plots and other analyzed results

```
cd figures_tables  
Rscript ../scripts/R/old_metaplots.R
Rscript ../scripts/R/new_metaplots.R 
Rscript ../scripts/R/multiplot_figures.R 
```

**NOTE #2:** These scripts were written for use on the MSU HPCC. To run them on your computer or a different environment, you will need to change the header of each shell script (python and R scripts do not require changing) to something appropriate for your system. In each shell script, you will also need to either modify or delete these lines to a path appropriate for your system:
```
export PATH="$HOME/miniconda3/envs/Bnapus-polyploidy/bin:$PATH" 
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/Bnapus-polyploidy/lib:$LD_LIBRARY_PATH"
```
Specifically modify **$HOME/miniconda3** to correspond to where your conda installation is.

You may also need to delete this line, which is meant for submitted jobs on the cluster.

```
cd $PBS_O_WORKDIR
```

**NOTE #3:** The file "scripts/py/functions.py" contains a number of various python functions that I have written that can be used to analyze data from methylpy. 