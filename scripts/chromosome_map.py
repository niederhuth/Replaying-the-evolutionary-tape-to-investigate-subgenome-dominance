import os
import sys

functionsfile = '../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions
allc='allc_'+sys.argv[1]+'.tsv.gz'
genome_file='../ref/combined/combined.genome'
filter_chr=['37_Plastid','37_Mitochondria','ChrL','ChrC','ChrM','Pt','Mt']
context=['CG','CHG','CHH','CAA','CAT','CAC','CAG','CTA','CTT','CTC',
         'CTG','CCA','CCT','CCC','CCG','CGA','CGT','CGC','CGG']

functions.genome_window_methylation(allc,genome_file,output='results/genome_window_methylation.tsv',
	window_size=100000,stepsize=50000,cutoff=0,filter_chr=filter_chr,output_mC_counts=True)

