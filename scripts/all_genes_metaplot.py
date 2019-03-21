import os
import sys
import pandas as pd

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
allc='allc_'+sys.argv[1]+'.tsv.gz'
annotations='../../ref/annotations/combined.gff'
genome_file='../../ref/annotations/combined.genome'
filter_chr=['ChrL','37_Plastid','37_Mitochondria']
mc_type=['CG','CHG','CHH']
window_number=60
updown_stream=2000
cutoff=0
first_feature='gene'
second_feature='CDS'
output='results/uncorrected_all_genes_metaplot.txt'

#get chromosome list
chrs = list(pd.read_table(genome_file,header=None,usecols=[0],dtype='str')[0])
chrs = list(set(chrs).difference(filter_chr))

#get gene metaplot data
functions.gene_metaplot(allc,annotations,genome_file,output=output,
	mc_type=mc_type,window_number=window_number,updown_stream=updown_stream,
	cutoff=cutoff,first_feature=first_feature,second_feature=second_feature,
	chrs=chrs,remove_tmp=False)
