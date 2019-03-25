import os
import sys
import pandas as pd

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
allc='uncorrected_all_genes_filered_allc.tmp'
annotations='../../ref/annotations/R500.gff'
genome_file='../../ref/combined/combined.genome'
filter_chr=['ChrL','37_Plastid','37_Mitochondria']
mc_type=['CG','CHG','CHH']
window_number=60
updown_stream=2000
cutoff=0
first_feature='gene'
second_feature='CDS'
filtered_data_output='uncorrected_R500_genes_filtered_allc.tmp'
output='results/uncorrected_R500_genes_metaplot.txt'

#get chromosome list
chrs = list(pd.read_table(genome_file,header=None,usecols=[0],dtype='str')[0])
chrs = [a for a in chrs if a[1] == "r"]

#get gene metaplot data
functions.gene_metaplot(allc,annotations,genome_file,output=output,
	mc_type=mc_type,window_number=window_number,updown_stream=updown_stream,
	cutoff=cutoff,first_feature=first_feature,second_feature=second_feature,
	chrs=chrs,filtered_data_output=filtered_data_output,remove_tmp=True)
