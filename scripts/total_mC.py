import os
import sys

functionsfile = '../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions
allc='allc_'+sys.argv[1]+'.tsv.gz'
genome_file1='../ref/TO1000/TO1000Ref.fa.fai'
genome_file2='../ref/IMB218/IMB218-CallR500.fasta.fai'
filter_chr=['37_Plastid','37_Mitochondria','ChrL','ChrC','ChrM','Pt','Mt']
context=['CG','CHG','CHH','CAA','CAT','CAC','CAG','CTA','CTT','CTC',
         'CTG','CCA','CCT','CCC','CCG','CGA','CGT','CGC','CGG']

functions.weighted_mC(allc, output="results/Total_methylation.tsv", cutoff=0)
functions.weighted_mC(allc, output="results/TO1000_total_methylation.tsv", cutoff=0, genome=genome_file1)
functions.weighted_mC(allc, output="results/IMB218_total_methylation.tsv", cutoff=0, genome=genome_file2)
