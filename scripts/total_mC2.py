import os
import sys
import pandas as pd

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
allc='corrected_allc.tsv'
genome_file='../../ref/combined/combined.genome'
genome_file2='../../ref/R500/R500.genome'
genome_file3='../../ref/TO1000/TO1000.genome'
filter_chr=['ChrL','37_Plastid','37_Mitochondria']
mc_type=['CG','CHG','CHH']
cutoff=0
output='results/corrected_total_weighted_methylation.txt'
output2='results/corrected_R500_weighted_methylation.txt'
output3='results/corrected_TO1000_weighted_methylation.txt'

#get chromosome list
chrs = list(pd.read_csv(genome_file,header=None,usecols=[0],dtype='str',sep="\t")[0])
chrs = list(set(chrs).difference(filter_chr))

#get total weighted mC
functions.total_weighted_mC(allc,output=output,mc_type=mc_type,cutoff=cutoff,chrs=chrs)

#get IMB218 chromosome list
chrs = list(pd.read_csv(genome_file2,header=None,usecols=[0],dtype='str',sep="\t")[0])
chrs = list(set(chrs).difference(filter_chr))

#get IMB218 total weighted mC
functions.total_weighted_mC(allc,output=output2,mc_type=mc_type,cutoff=cutoff,chrs=chrs)

#get TO1000 chromosome list
chrs = list(pd.read_csv(genome_file3,header=None,usecols=[0],dtype='str',sep="\t")[0])
chrs = list(set(chrs).difference(filter_chr))

#get TO1000 weighted mC
functions.total_weighted_mC(allc,output=output3,mc_type=mc_type,cutoff=cutoff,chrs=chrs)
