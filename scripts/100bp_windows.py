import os
import sys
import pandas as pd

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
allc='allc_'+sys.argv[1]+'.tsv.gz'
genome_file='../../ref/combined/combined.genome'
filter_chr=['ChrL','37_Plastid','37_Mitochondria']
mc_type=['CG','CHG','CHH']
window_size=100
stepsize=100
cutoff=0
output='results/100bp_windows.txt'

#get chromosome list
chrs = list(pd.read_csv(genome_file,header=None,usecols=[0],dtype='str',sep="\t")[0])
chrs = list(set(chrs).difference(filter_chr))

#get gene metaplot data
functions.genome_window_methylation(allc,genome_file,output=output,mc_type=mc_type,
	window_size=window_size,stepsize=stepsize,cutoff=cutoff,chrs=chrs)
