import os
import sys

functionsfile = '../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions
allc='allc_'+sys.argv[1]+'.tsv.gz'
genome_file='../ref/combined/combined.genome'
features='../ref/combined/combined.gff.gz'
filter_chr=['37_Plastid','37_Mitochondria','ChrL','ChrC','ChrM','Pt','Mt']
context=['CG','CHG','CHH','CAA','CAT','CAC','CAG','CTA','CTT','CTC',
         'CTG','CCA','CCT','CCC','CCG','CGA','CGT','CGC','CGG']

if os.path.exists(features):
	print('Getting gene methylation')
	functions.map2features(allc,features,genome_file,updown_stream=0,
        first_feature='gene',second_feature='CDS',filter_chr=filter_chr)

	functions.feature_mC_levels('CDS_allc.tmp',features,
		output='results/gene_methylation_levels.tsv',
        cutoff=0,filter_features='gene',filter_chr=filter_chr)
else:
	print('No annotations found')
