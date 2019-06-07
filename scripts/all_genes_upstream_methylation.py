import os
import sys
import pandas as pd
import pybedtools as pbt
from subprocess import call

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
allc='allc_'+sys.argv[1]+'.tsv.gz'
annotations='../../ref/annotations/combined.gff'
genome_file='../../ref/combined/combined.genome'
filter_chr=['ChrL','37_Plastid','37_Mitochondria']
mc_type=['CG','CHG','CHH']
updown_stream=2000
cutoff=0
first_feature='gene'
flank_bed='bed.tmp'
filtered_output="all_genes_filtered_allc.tmp"
output='results/all_genes_upstream_methylation.txt'

#get chromosome list
chrs = list(pd.read_csv(genome_file,header=None,usecols=[0],dtype='str',sep="\t")[0])
chrs = list(set(chrs).difference(filter_chr))

#extract upstream regions
print('Making gff file of flanking regions')
bed = pbt.BedTool(annotations)
pbt.bedtool.BedTool.flank(bed,g=genome_file,l=updown_stream,r=0,s=True).saveas('bed.tmp')
#correct sites where the start is greater than the end
command="awk -v FS='\t' -v OFS='\t' '{$5=($4>$5&&$4==1?$4:$5)}; {$4=($4>$5&&$4!=1?$5:$4)} 1' bed.tmp > tmp; mv tmp bed.tmp"
call(command, shell=True)

#get promoter methylation data
print('Getting upstream methylation data')
functions.feature_methylation(filtered_output,flank_bed,genome_file,output=output,
	mc_type=mc_type,updown_stream=0,
	feature=first_feature,cutoff=cutoff,chrs=chrs)

#remove tmp files
os.remove('bed.tmp')
