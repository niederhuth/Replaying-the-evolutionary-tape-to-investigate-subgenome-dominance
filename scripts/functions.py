import sys
import pandas as pd
import pybedtools as pbt
from os import remove
from numpy import float64
from subprocess import call
from io import TextIOWrapper
from itertools import product
from gzip import open as gzopen

#split large files into temporary files of smaller size
def split_file(input,line_number=10000000):
	#set starting count values
	a = 1
	b = 0
	c = line_number
	#open first output file and add to list
	files = ['split' + str(a) + '.tmp']
	out=open('split' + str(a) + '.tmp','w')
	#open input file
	with TextIOWrapper(gzopen(input,'rb')) as e:
		#iterate over each line, adding an index
		for index, line in enumerate(e):
			#test if index fits between upper and lower limits and write to file
			if index <= c:
				if index > b:
					out.write(str(line))
			else:
				#close last ouput
				out.close()
				#reset count values
				a += 1
				b = c
				c += line_number
				#open new output and add to list
				filesappend('split' + str(a) + '.tmp')
				out=open('split' + str(a) + '.tmp','w')
				#output line
				out.write(str(line))
		#close last ouptut
		out.close()
	#return number of temporary files for use in other functions
	return(files)

#function for filtering annotation files based on feature (gene, exon, mRNA, etc)
def feature_filter(x,feature):
	if feature:
		return x[2] == feature
	else:
		return x

#function for filtering annotation files based on strand
def strand_filter(x,strand):
	return x.strand == strand

#function for filtering annotation files based on chromosome
def chr_filter(x,chr):
	return x.chrom in chr

#interpret sequence context, taken from methylpy utilities
def expand_nucleotide_code(mc_type=['C']):
	iub_dict = {'N':['A','C','G','T','N'],'H':['A','C','T','H'],'D':['A','G','T','D'],'B':['C','G','T','B'],'A':['A','C','G','A'],'R':['A','G','R'],'Y':['C','T','Y'],'K':['G','T','K'],'M':['A','C','M'],'S':['G','C','S'],'W':['A','T','W'],'C':['C'],'G':['G'],'T':['T'],'A':['A']}
	mc_class = list(mc_type) # copy
	if 'C' in mc_type:
		mc_class.extend(['CGN', 'CHG', 'CHH','CNN'])
	elif 'CG' in mc_type:
		mc_class.extend(['CGN'])
	mc_class_final = []
	for motif in mc_class:
		mc_class_final.extend([''.join(i) for i in product(*[iub_dict[nuc] for nuc in motif])])
	return(set(mc_class_final))

#Read allc file and convert to bedfile
def allc2bed(allc,return_bed=True):
	#get first line
	if allc.endswith('gz'):
		header = gzopen(allc).readline().rstrip()
	else:
		header = open(allc).readline().rstrip()
	#check if first line contains header
	if header == 'chr\tpos\tstrand\tmc_class\tmc_count\ttotal\tmethylated':
		#read in allc file to pandas dataframe
		a = pd.read_csv(allc,dtype={'chr':str,'pos':int,'strand':str,'mc_class':str,'mc_count':int,'total':int,'methylated':int},sep="\t")
	else:
		#read in allc file to pandas dataframe, add header if does not have one
		a = pd.read_csv(allc,names=['chr','pos','strand','mc_class','mc_count','total','methylated'],dtype={'chr':str,'pos':int,'strand':str,'mc_class':str,'mc_count':int,'total':int,'methylated':int},sep="\t")
	#if return_bed = True, convert to bedfile
	if return_bed is True:
		#add new columns
		a['pos2'] = a['pos']
		a['name'] = a.index
		a['score'] = '.'
		#reorder columns
		a = a[['chr','pos','pos2','name','score','strand','mc_class','mc_count','total','methylated']]
		#create pybedtools object
		a = pbt.BedTool.from_dataframe(a)
	#return a
	return a

#Collect mC data for a context
def get_mC_data(a,mc_type='C',cutoff=0):
	#expand nucleotide list for a given context
	b = expand_nucleotide_code(mc_type=[mc_type])
	d1 = d2 = d3 = d4 = 0
	#iterate over each line
	for c in a.itertuples():
		#check if line is correct context
		if c[4] in b:
			#check if meets minimum cutoff for read depth
			if int(c[6]) >= int(cutoff):
				#add up number of sites
				d1 = d1 + 1
				#add up number of sites called methylated by methylpy
				d2 = d2 + int(c[7])
				#add up total reads covering a site
				d3 = d3 + int(c[6])
				#add up total methylated reads covering a site
				d4 = d4 + int(c[5])
	#create list
	e = [mc_type,d1,d2,d3,d4]
	#return that list
	return e

#Collect total methylation data for genome or sets of chromosomes
def total_weighted_mC(allc,output=(),mc_type=['CG','CHG','CHH'],cutoff=0,chrs=[]):
	#read allc file
	a = allc2bed(allc,return_bed=False)
	#filter chromosome sequences
	if chrs:
		a = a[a.chr.isin(chrs)]
	#create data frame
	columns=['Context','Total_sites','Methylated_sites','Total_reads','Methylated_reads','Weighted_mC']
	b = pd.DataFrame(columns=columns)
	#iterate over each mC type and run get_mC_data
	for c in mc_type:
		d = get_mC_data(a,mc_type=c,cutoff=cutoff)
		#calculate weighted methylation
		d += [(float64(d[4])/float64(d[3]))]
		#add to data frame
		b = b.append(pd.DataFrame([d],columns=columns), ignore_index=True)
	#output results
	if output:
		b.to_csv(output, sep='\t', index=False)
	else:
		return b

#for calculating methylation levels in windows across the genome
def genome_window_methylation(allc,genome_file,output=(),mc_type=['CG','CHG','CHH'],window_size=100000,stepsize=50000,cutoff=0,chrs=[]):
	#read in allc file
	print("Reading allc file")
	a = allc2bed(allc)
	#create output data frame
	print("Creating output dataframe")
	c = ['Chr','Window']
	columns=['Total_sites','Methylated_sites','Total_reads','Methylated_reads','Weighted_mC']
	for d in mc_type:
		for e in columns:
			c = c + [d + '_' + e]
	b = pd.DataFrame(columns=c)
	#make windows
	print("Making genome windows")
	w_bed = pbt.bedtool.BedTool.window_maker(pbt.BedTool(genome_file).filter(chr_filter,chrs),g=genome_file,w=window_size,s=stepsize,i='srcwinnum')
	#intersect bedfiles with pybedtools
	print("Mapping DNA methylation data to windows")
	mapping = pbt.bedtool.BedTool.intersect(a,w_bed,wa=True,wb=True)
	del(w_bed,a)
	#convert to pandas dataframe
	print("Converting to pandas dataframe")
	m = pd.read_csv(mapping.fn,header=None,usecols=[13,6,7,8,9],sep="\t")
	del(mapping)
	#split srcwinnum
	print("Formatting names")
	f = m[13].str.split('_', n = 1, expand = True)
	#make new columns from srcwinnum
	m['Chr'] = f[0]
	m['Window'] = f[1]
	del(f)
	#reorder data frame
	m = m[['Chr','Window',13,6,7,8,9]]
	#iterate over each chromosome
	print("Calculating window methylation data")
	for g in chrs:
		#get windows for that specific chromosome
		windows = list(m[m['Chr'].isin([str(g)])]['Window'].drop_duplicates())
		#iterate over each window in each chromosome
		for h in windows:
			#filter for rows matching specific chr & window number
			i = m[m['Chr'].isin([str(g)]) & m['Window'].isin([str(h)])]
			#make list for methylation data
			j = [g,h]
			#iterate over each mC type and run get_mC_data
			for k in mc_type:
				l = get_mC_data(i,mc_type=k,cutoff=cutoff)
				#delete first line of list
				del(l[0])
				#Calculate weighted methylation and add this to list of data for other mc_types
				j += l + [(float64(l[3])/float64(l[2]))]
			#append the results for that window to the dataframe
			b = b.append(pd.DataFrame([j],columns=c),ignore_index=True)
	#output results
	print("Outputting results")
	if output:
		b.to_csv(output, sep='\t', index=False)
	else:
		return b

#create filtered allc of sites mapping to annotations
def allc_annotation_filter(allc,annotations,genome_file,output=(),updown_stream=2000,first_feature=(),second_feature=(),chrs=[]):
	#read in annotations and filter by first feature, typically a something like 'gene'
	#this is used solely to accurately create flanking regions
	print("Reading annotations")
	bed = pbt.BedTool(annotations).filter(feature_filter,first_feature).filter(chr_filter,chrs)
	#create bedfile of flanking regions (if specified)
	print("Getting flanking regions")
	flank_bed = pbt.bedtool.BedTool.flank(bed,g=genome_file,l=updown_stream,r=updown_stream,s=True).saveas('f_bed.tmp')
	#correct sites where the start is greater than the end
	command="awk -v OFS='\t' '{$5=($4>$5&&$4==1?$4:$5)}; {$4=($4>$5&&$4!=1?$5:$4)} 1' f_bed.tmp > tmp; mv tmp f_bed.tmp"
	call(command, shell=True)
	#read in annotations and filter by second annotation, typically a something like coding sequences 'CDS'
	#this is the annotation used to first filter the data
	print("Filtering second feature annotations")
	cds_bed = pbt.BedTool(annotations).filter(feature_filter,second_feature).filter(chr_filter,chrs).saveas('c_bed.tmp')
	#combine flanking regions with second feature
	print("Combining second feature and flanking region bed files")
	combined_bed = cds_bed.cat(flank_bed, postmerge=False)
	#read in allc file and map to annotations
	print("Reading allc file")
	a = allc2bed(allc)
	print("Mapping sites to annotations")
	mapping = pbt.bedtool.BedTool.intersect(a,combined_bed,wa=True)
	print("Converting mapped sites to table")
	m = pd.read_csv(mapping.fn, header=None, usecols = [0,1,5,6,7,8,9],sep="\t")
	#create new filtered allc file of sites mapping to regions
	print("Reformat data and drop duplicate sites")
	m.columns = ['chr','pos','strand','mc_class','mc_count','total','methylated']
	m = m.drop_duplicates().sort_values(['chr','pos'],ascending=[True,True])
	#output results	#output results
	if output:
		print("Outputing results")
		m.to_csv(output, sep='\t', index=False)
	else:
		return m
	#remove temporary files created
	print("Removing temporary files")
	tmp=['f_bed.tmp','c_bed.tmp']
	for b in tmp:
		remove(b)

#output methylation data for making metaplots of features (e.g. genes, CDS, transposons), will not filter out data from introns, etc...use gene_metaplot for that
def metaplot(allc,annotations,genome_file,output=(),mc_type=['CG','CHG','CHH'],window_number=60,updown_stream=2000,feature=(),cutoff=0,chrs=[]):
	#read in allc file
	print("Reading allc file")
	a = allc2bed(allc)
	#create output data frame
	print("Create output dataframe")
	c = ['Window']
	columns=['Weighted_mC']
	for d in mc_type:
		for e in columns:
			c = c + [d + '_' + e]
	b = pd.DataFrame(columns=c)
	#read annotation file and filter on feature
	print("Reading annotations")
	f_bed = pbt.BedTool(annotations).filter(feature_filter,feature).filter(chr_filter,chrs).saveas('f_bed.tmp')
	#if updown_stream set to 0, only region to look at is the feature itself, e.g. f_bed
	if updown_stream == 0:
		print("No flanking regions specified. Analyze annotations only")
		regions=[f_bed]
	#if updown_stream specified, create bed files for upstream regions (u_bed) and down stream regions (d_bed)
	else:
		print("Get flanking regions")
		u_bed = pbt.bedtool.BedTool.flank(f_bed,g=genome_file,l=updown_stream,r=0,s=True).saveas('u_bed.tmp')
		command="awk -v OFS='\t' '{$5=($4>$5&&$4==1?$4:$5)}; {$4=($4>$5&&$4!=1?$5:$4)} 1' u_bed.tmp > tmp; mv tmp u_bed.tmp"
		call(command, shell=True)
		d_bed = pbt.bedtool.BedTool.flank(f_bed,g=genome_file,l=0,r=updown_stream,s=True).saveas('d_bed.tmp')
		command="awk -v OFS='\t' '{$5=($4>$5&&$4==1?$4:$5)}; {$4=($4>$5&&$4!=1?$5:$4)} 1' d_bed.tmp > tmp; mv tmp d_bed.tmp"
		call(command, shell=True)
		regions=[u_bed,f_bed,d_bed]
	#set window number to 1
	window = 1
	#iterate over each region and collect methylation data
	print("Mapping sites to annotations and collecting methylation data")
	for f in regions:
		#filter bed files based on strand
		p_bed = f.filter(strand_filter,strand='+').saveas('p_bed.tmp')
		n_bed = f.filter(strand_filter,strand='-').saveas('n_bed.tmp')
		#make windows
		pw_bed = pbt.bedtool.BedTool.window_maker(p_bed,b=p_bed,n=int(window_number/3),i='srcwinnum')
		nw_bed = pbt.bedtool.BedTool.window_maker(n_bed,b=n_bed,n=int(window_number/3),i='srcwinnum',reverse=True)
		#combine windows
		w_bed = pw_bed.cat(nw_bed, postmerge=False)
		#intersect bedfiles with pybedtools
		mapping = pbt.bedtool.BedTool.intersect(a,w_bed,wa=True,wb=True)
		del(w_bed,pw_bed,nw_bed)
		#convert to pandas dataframe
		m = pd.read_csv(mapping.fn,header=None,usecols=[13,6,7,8,9],sep="\t")
		del(mapping)
		#split srcwinnum
		g = m[13].str.split('_', n = 1, expand = True)
		#make new columns from srcwinnum
		m['Name'] = g[0]
		m['Window'] = g[1]
		del(g)
		#reorder data frame
		m = m[['Name','Window',13,6,7,8,9]]
		#iterate list of window numbers
		for h in list(range(1,int(window_number/3)+1)):
			#filter for rows matching specific window number
			i = m[m['Window'].isin([str(h)])]
			#make list for methylation data
			j = [window]
			#iterate over each mC type and run get_mC_data
			for k in mc_type:
				l = get_mC_data(i,mc_type=k,cutoff=cutoff)
				#Calculate weighted methylation and add this to list of data for other mc_types
				j += [(float64(l[4])/float64(l[3]))]
			#append the results for that window to the dataframe
			b = b.append(pd.DataFrame([j],columns=c),ignore_index=True)
			#count windows
			window += 1
	#output results
	if output:
		print("Outputing results")
		b.to_csv(output, sep='\t', index=False)
	else:
		return b
	#remove temporary files created
	print("Removing temporary files")
	tmp=['f_bed.tmp','u_bed.tmp','d_bed.tmp','p_bed.tmp','n_bed.tmp']
	for n in tmp:
		remove(n)

#output methylation data for making metaplots of features (e.g. genes, CDS, transposons), for use when need to first filter data from another feature, such as intron sequences
def gene_metaplot(allc,annotations,genome_file,output=(),mc_type=['CG','CHG','CHH'],window_number=60,updown_stream=2000,cutoff=0,first_feature=(),second_feature=(),chrs=[],filtered_data_output='annotation_filtered_allc.tmp',remove_tmp=True):
	#prefilter allc file based on annotations
	print("Filtering data for sites in annotations")
	allc_annotation_filter(allc,annotations,genome_file,output=filtered_data_output,updown_stream=updown_stream,first_feature=first_feature,second_feature=second_feature,chrs=chrs)
	#collect methylation data
	print("Collecting metaplot data")
	metaplot(filtered_data_output,annotations,genome_file,output=output,mc_type=mc_type,window_number=window_number,updown_stream=updown_stream,feature=first_feature,cutoff=0,chrs=chrs)
	#remove annotation filtered allc file, if set to false, this will be kept
	if remove_tmp:
		remove(filtered_data_output)

#Calculate methylation levels for features
def feature_methylation(allc,annotations,genome_file,output=(),mc_type=['CG','CHG','CHH'],updown_stream=0,feature=(),cutoff=0,chrs=[]):
	#read in allc file
	a = allc2bed(allc)
	#create output data frame
	c = ['Feature']
	columns=['Total_C','Methylated_mC','Total_Reads','Methylated_Reads','Weighted_mC']
	for d in mc_type:
		for e in columns:
			c = c + [d + '_' + e]
	b = pd.DataFrame(columns=c)
	#read annotation file and filter on feature
	f_bed = pbt.BedTool(annotations).filter(feature_filter,feature).filter(chr_filter,chrs).saveas('f_bed.tmp')
	#if updown_stream set to 0, only region to look at is the feature itself, e.g. f_bed
	if updown_stream == 0:
		regions=['f_bed.tmp']
	#if updown_stream specified, create bed files for upstream regions (u_bed) and down stream regions (d_bed)
	else:
		u_bed = pbt.bedtool.BedTool.flank(f_bed,g=genome_file,l=updown_stream,r=0,s=True).saveas('u_bed.tmp')
		command="awk -v OFS='\t' '{$5=($4>$5&&$4==1?$4:$5)}; {$4=($4>$5&&$4!=1?$5:$4)} 1' u_bed.tmp > tmp; mv tmp u_bed.tmp"
		call(command, shell=True)
		d_bed = pbt.bedtool.BedTool.flank(f_bed,g=genome_file,l=0,r=updown_stream,s=True).saveas('d_bed.tmp')
		command="awk -v OFS='\t' '{$5=($4>$5&&$4==1?$4:$5)}; {$4=($4>$5&&$4!=1?$5:$4)} 1' d_bed.tmp > tmp; mv tmp d_bed.tmp"
		regions=['u_bed.tmp','d_bed.tmp']
	#iterate over each region and collect methylation data
	for f in regions:
		#intersect bedfiles with pybedtools
		mapping = pbt.bedtool.BedTool.intersect(a,f,wa=True,wb=True)
		del(a)
		#convert to pandas dataframe
		m = pd.read_csv(mapping.fn,header=None,usecols=[18,6,7,8,9],sep="\t")
		del(mapping)
		#split column 18
		g = m[18].str.split(';', n = 1, expand = True)
		g = g[0].str.split('=', n = 1, expand = True)
		#make new columns from column 18
		m['Name'] = g[1]
		m['Window'] = '.'
		#reorder m
		m = m[['Name','Window',18,6,7,8,9]]
		#get gene names
		g = m['Name'].drop_duplicates()
		#iterate list of gene names
		for h in list(g):
			#filter for rows matching specific gene
			i = m[m['Name'].isin([str(h)])]
			#make list for methylation data
			j = [h]
			#iterate over each mC type and run get_mC_data
			for k in mc_type:
				#check if i is empty
				if not i:
					#if empty, add 0 column
					j += ['NA','NA','NA','NA','NA']
				#else if not empty
				else:
					l = get_mC_data(i,mc_type=k,cutoff=cutoff)
					#Calculate weighted methylation and add this to list of data for other mc_types
					j += [l[1],l[2],l[3],l[4],(float64(l[4])/float64(l[3]))]
			#append the results for that window to the dataframe
			b = b.append(pd.DataFrame([j],columns=c),ignore_index=True)
	#output results
	if output:
		b.to_csv(output, sep='\t', index=False)
	else:
		return b
	#remove temporary files created
	for n in regions:
		remove(n)
