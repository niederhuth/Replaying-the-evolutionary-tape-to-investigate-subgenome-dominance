import os
import sys

functionsfile = '../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions
allc='allc_'+sys.argv[1]+'.tsv.gz'
genome_file='../ref/combined/combined.genome'
features='../ref/combined/'+sys.argv[1]+'.gff.gz'
filter_chr=['37_Plastid','37_Mitochondria','ChrL','ChrC','ChrM','Pt','Mt']
context=['CG','CHG','CHH','CAA','CAT','CAC','CAG','CTA','CTT','CTC',
         'CTG','CCA','CCT','CCC','CCG','CGA','CGT','CGC','CGG']

if os.path.exists(genome_file):
	print('Getting window methylation')
	functions.genome_window_methylation(allc,genome_file,
		output='results/genome_window_methylation.tsv',window_size=100,cutoff=0,
		filter_chr=filter_chr,output_mC_counts=True)

	#print('Getting gene methylation')
	#functions.map2features(allc,features,genome_file,updown_stream=0,
        #	first_feature='gene',second_feature='CDS',filter_chr=filter_chr)
	
	#functions.feature_mC_levels('CDS_allc.tmp',features,
	#	output='results/gene_methylation_levels.tsv',
	#	cutoff=0,filter_features='gene',filter_chr=filter_chr)

	#functions.gene_binom_test('results/gene_methylation_levels.tsv',cutoff=10,
	#	calc_baseline=False,mCG=0.17390540071880004,mCHG=0.04645834951454248,mCHH=0.004470085236457535,
	#	output='results/binom_test.tsv')
	#functions.classify_genes('results/binom_test.tsv',output='results/classified_genes.tsv',min_sites=0,qvalue=0.05)

	#os.remove('CDS_allc.tmp')	
	#os.remove('f_tmp')
	#os.remove('c_tmp')
else:
	print('No annotations found')
