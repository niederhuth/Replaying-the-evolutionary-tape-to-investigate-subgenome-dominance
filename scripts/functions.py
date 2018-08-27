import os
import sys
import io
import pybedtools as pbt
import pandas as pd
import itertools
import numpy as np
import gffutils
import gzip
from scipy import stats
from scipy.stats.stats import pearsonr
from collections import Counter
from Bio import SeqIO
from math import ceil

#For splitting a dataframe so that columns of equal value do not get put in seperate files
def split_df_on_column(m,size=10000000,column=0):
	f=ceil(len(m)/size)
	w=0
	x=size
	y=size
	for i in range(1,f):
		z=True
		while z:
			if m.iloc[x][column] == m.iloc[y][column]:
				y=y+1
			else:
				z=False
		m.iloc[w:y-1].to_csv('tmp'+str(i), sep='\t', index=False)
		w=y
		x=y+size
		y=y+size
	m.iloc[w:].to_csv('tmp'+str(f), sep='\t', index=False)
	return(f)

#For splitting large files
def split_large_file(input,lines=10000000):
    with gzip.open(input,'rb') as a:
        l=0
        for line in a:
            l=l+1
    f=ceil(l/lines)
    x=0
    y=lines
    for i in range(1,f+1):
        out=open('tmp'+str(i),'w')
        with io.TextIOWrapper(gzip.open(input,'rb')) as a:
            for index, line in enumerate(a):
                if index <= y:
                    if index > x:
                        out.write(str(line))
                else:
                    break
        out.close()
        x=y+1
        y=y+lines+1
    return(f)

#interpret sequence context, taken from methylpy.utils
def expand_nucleotide_code(mc_type=["C"]):
    iub_dict = {"N":["A","C","G","T"],"H":["A","C","T"],"C":["C"],"G":["G"],"T":["T"],"A":["A"]}
    for type in mc_type[:]:
        type += "N" * (3 - len(type))
        mc_type.extend(["".join(i) for i in itertools.product(*[iub_dict[nuc] for nuc in type])])
    if "C" in mc_type:
        mc_type.extend(["CG", "CHG", "CHH","CNN"])
    if "CG" in mc_type:
        mc_type.extend(["CGN"])
    return mc_type

#filter allc file based on sequence context
def filter_context(allc,context=["C"]):
	if allc.endswith('gz'):
		header = gzip.open(allc).readline().rstrip()
	else:
		header = open(allc).readline().rstrip()
	if header == 'chr\tpos\tstrand\tmc_class\tmc_count\ttotal\tmethylated':
		a = pd.read_table(allc,dtype={'chr':str,'pos':int,'strand':str,'mc_class':str,
			'mc_count':int,'total':int,'methylated':int})
	else:
		a = pd.read_table(allc,names=['chr','pos','strand','mc_class','mc_count',
			'total','methylated'],dtype={'chr':str,'pos':int,'strand':str,
			'mc_class':str,'mc_count':int,'total':int,'methylated':int})
	a = a[a.mc_class.isin(expand_nucleotide_code(context))]
	return a

#converts and allc file to bed format for use with bedtools
def allc2bed(allc,context=["C"],bed=True):
    a = filter_context(allc,context)
    a['pos2'] = a.pos
    a['name'] = a.index
    a['score'] = "."
    a = a[['chr','pos','pos2','name','score','strand','mc_class','mc_count','total','methylated']]
    if bed is True:
        a = pbt.BedTool.from_dataframe(a)
    return a

#simple function for filtering gff files based on feature (gene, exon, mRNA, etc)
def feat_filter(x,feature):
    if feature:
        return x[2] == feature
    else:
        return x

#simple function for filtering gff files based on strand
def strand_filter(x,strand):
    return x.strand == strand

#simple function for filtering gff files based on chromosome
def chr_filter(x,chr):
    return x.chrom not in chr

#make an intron gff file
def intron_gff(features,output=(),remove_db=True):
    db=gffutils.create_db(features, 'gff.db')
    gff_out=gffutils.gffwriter.GFFWriter(output, in_place=True)
    name="none"
    count=1
    for rec in db.create_introns():
    	if rec['Parent'] == name:
        	count=count+1
        	rec.attributes['ID'] = [str(rec['Parent'])[2:-2] + '.intron' + str(count)]
       		gff_out.write_rec(rec)
    	else:
        	name = rec['Parent']
        	count = 1
        	rec.attributes['ID'] = [str(rec['Parent'])[2:-2] + '.intron' + str(count)]
        	gff_out.write_rec(rec)
    gff_out.close()
    if remove_db:
        os.remove('gff.db')

#generic function for calculating methylation levels in windows
def window_methylation_levels(m,cutoff=0,nuc_bed=(),output_mC_counts=False):
    if nuc_bed:
        if output_mC_counts:
            a = pd.DataFrame(columns=['window','mCG_reads','CG_reads','mCG',
                                      'mCHG_reads','CHG_reads','mCHG','mCHH_reads',
                                       'CHH_reads','mCHH','coverage'])
        else:
             a = pd.DataFrame(columns=['window','mCG','mCHG','mCHH','coverage'])
    else:
        if output_mC_counts:
            a = pd.DataFrame(columns=['window','mCG_reads','CG_reads','mCG',
                                      'mCHG_reads','CHG_reads','mCHG','mCHH_reads',
                                       'CHH_reads','mCHH'])
        else:
             a = pd.DataFrame(columns=['window','mCG','mCHG','mCHH'])
    name = "none"
    C = CG = mCG = CHG = mCHG = CHH = mCHH = GC = 0
    if nuc_bed:
        nuc = pd.read_table(nuc_bed.fn, usecols = [3,7,8])
        m = pd.merge(m,nuc,left_on=13,right_on='4_usercol')
    for c in m.itertuples():
      if name == "none":
        name = c[5]
        if nuc_bed:
            GC = int(c[7]) + int(c[8])
        if int(c[3]) >= int(cutoff):
          C = C + 1
          if c[1].startswith("CN") or c[1].endswith("N"):
            continue
          elif c[1].startswith("CG"):
            CG = CG + int(c[3])
            mCG = mCG + int(c[2])
          elif c[1].endswith("G"):
            CHG = CHG + int(c[3])
            mCHG = mCHG + int(c[2])
          else:
            CHH = CHH + int(c[3])
            mCHH = mCHH + int(c[2])
      elif c[5] != name:
        if nuc_bed:
            if output_mC_counts:
                a = a.append({'window':str(name), 'mCG_reads':int(mCG), 'CG_reads':int(CG), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG_reads':int(mCHG), 'CHG_reads':int(CHG), 'mCHG':(np.float64(mCHG)/np.float64(CHG)),  'mCHH_reads':int(mCHH), 'CHH_reads':int(CHH), 'mCHH':(np.float64(mCHH)/np.float64(CHH)), 'coverage':(np.float64(C)/np.float(GC))}, ignore_index=True)
            else:
                a = a.append({'window':str(name), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG':(np.float64(mCHG)/np.float64(CHG)), 'mCHH':(np.float64(mCHH)/np.float64(CHH)), 'coverage':(np.float64(C)/np.float(GC))}, ignore_index=True)
        else:
            if output_mC_counts:
                a = a.append({'window':str(name), 'mCG_reads':int(mCG), 'CG_reads':int(CG), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG_reads':int(mCHG), 'CHG_reads':int(CHG), 'mCHG':(np.float64(mCHG)/np.float64(CHG)),  'mCHH_reads':int(mCHH), 'CHH_reads':int(CHH), 'mCHH':(np.float64(mCHH)/np.float64(CHH))}, ignore_index=True)
            else:
                a = a.append({'window':str(name), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG':(np.float64(mCHG)/np.float64(CHG)), 'mCHH':(np.float64(mCHH)/np.float64(CHH))}, ignore_index=True)
        name = c[5]
        if nuc_bed:
            GC = int(c[7]) + int(c[8])
        C = CG = mCG = CHG = mCHG = CHH = mCHH = 0
        if int(c[3]) >= int(cutoff):
            C = C + 1
            if c[1].startswith("CN") or c[1].endswith("N"):
              continue
            elif c[1].startswith("CG"):
              CG = CG + int(c[3])
              mCG = mCG + int(c[2])
            elif c[1].endswith("G"):
              CHG = CHG + int(c[3])
              mCHG = mCHG + int(c[2])
            else:
              CHH = CHH + int(c[3])
              mCHH = mCHH + int(c[2])
      elif c[5] == name:
        C = C + 1
        if int(c[3]) >= int(cutoff):
          if c[1].startswith("CN") or c[1].endswith("N"):
            continue
          elif c[1].startswith("CG"):
            CG = CG + int(c[3])
            mCG = mCG + int(c[2])
          elif c[1].endswith("G"):
            CHG = CHG + int(c[3])
            mCHG = mCHG + int(c[2])
          else:
            CHH = CHH + int(c[3])
            mCHH = mCHH + int(c[2])
    if nuc_bed:
        if output_mC_counts:
            a = a.append({'window':str(name), 'mCG_reads':int(mCG), 'CG_reads':int(CG), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG_reads':int(mCHG), 'CHG_reads':int(CHG), 'mCHG':(np.float64(mCHG)/np.float64(CHG)), 'mCHH_reads':int(mCHH), 'CHH_reads':int(CHH), 'mCHH':(np.float64(mCHH)/np.float64(CHH)), 'coverage':(np.float64(C)/np.float(GC))}, ignore_index=True)
        else:
            a = a.append({'window':str(name), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG':(np.float64(mCHG)/np.float64(CHG)), 'mCHH':(np.float64(mCHH)/np.float64(CHH)), 'coverage':(np.float64(C)/np.float(GC))}, ignore_index=True)
    else:
        if output_mC_counts:
            a = a.append({'window':str(name), 'mCG_reads':int(mCG), 'CG_reads':int(CG), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG_reads':int(mCHG), 'CHG_reads':int(CHG), 'mCHG':(np.float64(mCHG)/np.float64(CHG)), 'mCHH_reads':int(mCHH), 'CHH_reads':int(CHH), 'mCHH':(np.float64(mCHH)/np.float64(CHH))}, ignore_index=True)
        else:
            a = a.append({'window':str(name), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG':(np.float64(mCHG)/np.float64(CHG)), 'mCHH':(np.float64(mCHH)/np.float64(CHH))}, ignore_index=True)
    return a

#map methylation to features
def map2features(allc,features,genome_file,updown_stream=2000,first_feature=(),second_feature=(),filter_chr=[]):
	f = split_large_file(allc,lines=10000000)
	bed = pbt.BedTool(features).filter(feat_filter,first_feature).filter(chr_filter,filter_chr)
	flank_bed = pbt.bedtool.BedTool.flank(bed,g=genome_file,l=updown_stream,r=updown_stream,s=True).saveas('f_tmp')
	cds_bed = pbt.BedTool(features).filter(feat_filter,second_feature).filter(chr_filter,filter_chr).saveas('c_tmp')
	bed = cds_bed.cat(flank_bed, postmerge=False).saveas('b_tmp')
	tables=[]
	tables2=[]
	for i in range(1,f+1):
		tables.append('tmp'+str(i))
		mC_bed = allc2bed('tmp'+str(i))
		mapping = pbt.bedtool.BedTool.intersect(mC_bed,bed,wa=True).saveas('CDS'+str(i)+'.tmp')
		tables2.append('CDS'+str(i)+'.tmp')
		del(mapping)
	df_from_each_tmp_file = (pd.read_table(i,header=None) for i in tables2)
	m = pd.concat(df_from_each_tmp_file, ignore_index=True)
	m = m.drop_duplicates()
	m = m[[0,1,5,6,7,8,9]]
	m.columns = ['chr','pos','strand','mc_class','mc_count','total','methylated']
	m.to_csv('CDS_allc.tmp', sep='\t', index=False)
	del(m)
	for i in tables + tables2 + ['f_tmp','c_tmp', 'b_tmp']:
		os.remove(i)

#Get feature methylation data
def feature_mC_levels(allc,features,output=(),cutoff=0,filter_features=(),filter_chr=[]):
    bed = pbt.BedTool(features).filter(feat_filter,filter_features).filter(chr_filter,filter_chr)
    mC_bed = allc2bed(allc)
    mapping = pbt.bedtool.BedTool.intersect(mC_bed,bed,wa=True,wb=True)
    m = pd.read_table(mapping.fn, header=None, usecols = [18,6,7,8,9])
    m = m.sort_values(by = 18,ascending=True)
    a = pd.DataFrame(columns=['Gene','CG_sites','mCG_sites','CG_reads','mCG_reads',
                              'CG_methylation_level','CHG_sites','mCHG_sites','CHG_reads',
                              'mCHG_reads','CHG_methylation_level','CHH_sites','mCHH_sites',
                              'CHH_reads','mCHH_reads','CHH_methylation_level','C_sites',
                              'mC_sites','C_reads','mC_reads','C_methylation_level'])
    name = "none"
    rCG = mrCG = CG = mCG = rCHG = mrCHG = CHG = mCHG = rCHH = mrCHH = CHH = mCHH = 0
    for c in m.itertuples():
        if name == "none":
            name = c[5]
            if int(c[3]) >= int(cutoff):
                if c[1].startswith("CN"):
                    continue
                elif c[1].startswith("CG"):
                    rCG = rCG + int(c[3])
                    mrCG = mrCG + int(c[2])
                    CG = CG + 1
                    mCG = mCG + int(c[4])
                elif c[1].endswith("G"):
                    rCHG = rCHG + int(c[3])
                    mrCHG = mrCHG + int(c[2])
                    CHG = CHG + 1
                    mCHG = mCHG + int(c[4])
                else:
                    rCHH = rCHH + int(c[3])
                    mrCHH = mrCHH + int(c[2])
                    CHH = CHH + 1
                    mCHH = mCHH + int(c[4])
        elif c[5] != name:
            a = a.append({'Gene':name, 'CG_sites':CG, 'mCG_sites':mCG,
                          'CG_reads':rCG, 'mCG_reads':mrCG,
                          'CG_methylation_level':str(np.float64(mrCG)/np.float64(rCG)),
                          'CHG_sites':CHG, 'mCHG_sites':mCHG, 'CHG_reads':rCHG,
                          'mCHG_reads':mrCHG, 'CHG_methylation_level':str(np.float64(mrCHG)/np.float64(rCHG)),
                          'CHH_sites':CHH, 'mCHH_sites':mCHH, 'CHH_reads':rCHH, 'mCHH_reads':mrCHH,
                          'CHH_methylation_level':str(np.float64(mrCHH)/np.float64(rCHH)),
                          'C_sites':(CG+CHG+CHH), 'mC_sites':(mCG+mCHG+mCHH),
                          'C_reads':(rCG+rCHG+rCHH), 'mC_reads':(mrCG+mrCHG+mrCHH),
                          'C_methylation_level':str(np.float64(mrCG+mrCHG+mrCHH)/np.float64(rCG+rCHG+rCHH))
                          },ignore_index=True)
            name = c[5]
            rCG = mrCG = CG = mCG = rCHG = mrCHG = CHG = mCHG = rCHH = mrCHH = CHH = mCHH = 0
            if int(c[3]) >= int(cutoff):
                if c[1].startswith("CN"):
                    continue
                elif c[1].startswith("CG"):
                    rCG = rCG + int(c[3])
                    mrCG = mrCG + int(c[2])
                    CG = CG + 1
                    mCG = mCG + int(c[4])
                elif c[1].endswith("G"):
                    rCHG = rCHG + int(c[3])
                    mrCHG = mrCHG + int(c[2])
                    CHG = CHG + 1
                    mCHG = mCHG + int(c[4])
                else:
                    rCHH = rCHH + int(c[3])
                    mrCHH = mrCHH + int(c[2])
                    CHH = CHH + 1
                    mCHH = mCHH + int(c[4])
        elif c[5] == name:
            if int(c[3]) >= int(cutoff):
                if c[1].startswith("CN"):
                    continue
                elif c[1].startswith("CG"):
                    rCG = rCG + int(c[3])
                    mrCG = mrCG + int(c[2])
                    CG = CG + 1
                    mCG = mCG + int(c[4])
                elif c[1].endswith("G"):
                    rCHG = rCHG + int(c[3])
                    mrCHG = mrCHG + int(c[2])
                    CHG = CHG + 1
                    mCHG = mCHG + int(c[4])
                else:
                    rCHH = rCHH + int(c[3])
                    mrCHH = mrCHH + int(c[2])
                    CHH = CHH + 1
                    mCHH = mCHH + int(c[4])
    a = a.append({'Gene':name, 'CG_sites':CG, 'mCG_sites':mCG,
                  'CG_reads':rCG, 'mCG_reads':mrCG,
                  'CG_methylation_level':str(np.float64(mrCG)/np.float64(rCG)),
                  'CHG_sites':CHG, 'mCHG_sites':mCHG, 'CHG_reads':rCHG,
                  'mCHG_reads':mrCHG, 'CHG_methylation_level':str(np.float64(mrCHG)/np.float64(rCHG)),
                  'CHH_sites':CHH, 'mCHH_sites':mCHH, 'CHH_reads':rCHH, 'mCHH_reads':mrCHH,
                  'CHH_methylation_level':str(np.float64(mrCHH)/np.float64(rCHH)),
                  'C_sites':(CG+CHG+CHH), 'mC_sites':(mCG+mCHG+mCHH),
                  'C_reads':(rCG+rCHG+rCHH), 'mC_reads':(mrCG+mrCHG+mrCHH),
                  'C_methylation_level':str(np.float64(mrCG+mrCHG+mrCHH)/np.float64(rCG+rCHG+rCHH))
                  },ignore_index=True)
    if output:
        a.to_csv(output, sep='\t', index=False)
    else:
        return a

#
def FDR(a,column,new_col):
    b=a[a[column] != "NaN"][column].sort_values(ascending=False)
    lp=len(b)
    x=list(reversed(range(1,lp+1)))
    y=np.minimum.accumulate([lp/i*j for i,j in zip(x,b)])
    z=[i if i < 1.0 else 1.0 for i in y]
    c=pd.DataFrame(z,columns=[new_col])
    c.index=b.index
    a=pd.concat([a,c],axis=1)
    return a

#
def gene_binom_test(df,cutoff=10,calc_baseline=True,mCG=(),mCHG=(),mCHH=(),output=()):
    a=pd.read_table(df,sep="\t")
    if calc_baseline:
        mCG=pd.DataFrame.sum(a['mCG_sites'])/pd.DataFrame.sum(a['CG_sites'])
        mCHG=pd.DataFrame.sum(a['mCHG_sites'])/pd.DataFrame.sum(a['CHG_sites'])
        mCHH=pd.DataFrame.sum(a['mCHH_sites'])/pd.DataFrame.sum(a['CHH_sites'])
    elif not mCG or not mCHG or not mCHH:
        print('Use must specify baseline mCG, mCHG, and mCHH levels')
        return
    a['CG_pvalue']=stats.binom.sf(a.mCG_sites-1,a.CG_sites,mCG)
    a['CG_pvalue']=[i if s >= cutoff else "NaN" for i,s in zip(a['CG_pvalue'],a['CG_sites'])]
    a=FDR(a,column="CG_pvalue",new_col="CG_qvalue")
    a['CHG_pvalue']=stats.binom.sf(a.mCHG_sites-1,a.CHG_sites,mCHG)
    a['CHG_pvalue']=[i if s >= cutoff else "NaN" for i,s in zip(a['CHG_pvalue'],a['CHG_sites'])]
    a=FDR(a,column="CHG_pvalue",new_col="CHG_qvalue")
    a['CHH_pvalue']=stats.binom.sf(a.mCHH_sites-1,a.CHH_sites,mCHH)
    a['CHH_pvalue']=[i if s >= cutoff else "NaN" for i,s in zip(a['CHH_pvalue'],a['CHH_sites'])]
    a=FDR(a,column="CHH_pvalue",new_col="CHH_qvalue")
    if output:
        a.to_csv(output, sep='\t', index=False)
    else:
        return a

#
def classify(row,min_sites=0,qvalue=0.05):
    if row['CHH_qvalue'] <= qvalue and row['CHH_sites'] >= min_sites:
        return 'CHH'
    elif row['CHG_qvalue'] <= qvalue and row['CHG_sites'] >= min_sites:
        return 'CHG'
    elif row['CG_qvalue'] <= qvalue and row['CG_sites'] >= min_sites:
        return 'CG'
    elif (row['mCG_sites']+row['mCHG_sites']+row['mCHH_sites']) == 0 and \
         (row['CG_sites']+row['CHG_sites']+row['CHH_sites']) >= min_sites:
        return 'UM'
    else:
        return 'NA'

#
def classify_genes(df,output=(),min_sites=0,qvalue=0.05):
    a = pd.read_table(df,sep="\t")
    a['Classification'] = a.apply(classify,axis=1,min_sites=min_sites,qvalue=qvalue)
    if output:
        a.to_csv(output, sep='\t', index=False)
    else:
        return a

#
def feature_window_methylation(allc,features,output=(),window_size=100,filter_features=(),cutoff=0,filter_chr=[],output_mC_counts=True,ignore_strand=False):
    db=gffutils.create_db(features, 'gff.db')
    a=open('tmp.bed','w')
    for rec in db.all_features(featuretype=filter_features):
        a.write(str(rec[0])+'\t'+str(rec[3]-1)+'\t'+str(rec[4])+'\t'+
                str(rec.attributes['ID'])[2:-2]+'\t'+str(rec[5])+'\t'+
                str(rec[6])+'\n')
    a.close()
    bed = pbt.BedTool('tmp.bed').filter(chr_filter,filter_chr).saveas('b_tmp')
    if ignore_strand:
    	w_bed = pbt.bedtool.BedTool.window_maker(bed,b=bed,w=window_size,i='srcwinnum')
    else:
    	p_bed = bed.filter(strand_filter,strand='+').saveas('p_tmp')
    	n_bed = bed.filter(strand_filter,strand='-').saveas('n_tmp')
    	pw_bed = pbt.bedtool.BedTool.window_maker(p_bed,b=p_bed,w=window_size,i='srcwinnum')
    	nw_bed = pbt.bedtool.BedTool.window_maker(n_bed,b=n_bed,w=window_size,i='srcwinnum',reverse=True)
    	w_bed = pw_bed.cat(nw_bed, postmerge=False)
    mC_bed = allc2bed(allc)
    allc_mapping = pbt.bedtool.BedTool.intersect(mC_bed,w_bed,wa=True,wb=True)
    m = pd.read_table(allc_mapping.fn,header=None,usecols=[10,13,6,7,8])
    #m = m.sort_values(by = 13,ascending=True)
    b = window_methylation_levels(m,cutoff=cutoff,nuc_bed=(),output_mC_counts=True)
    if output:
        b.to_csv(output, sep='\t', index=False)
    else:
        return b
    for i in ['tmp.bed','gff.db','p_tmp','n_tmp','b_tmp']:
        os.remove(i)
    for i in ['w_bed','mC_bed','allc_mapping','m','b','pw_bed','nw_bed']:
        del(i)

#For calculating methylation levels in windows across the genome
def genome_window_methylation(allc,genome_file,output=(),window_size=100000,stepsize=50000,cutoff=0,filter_chr=[],output_mC_counts=True):
	f = split_large_file(allc,lines=50000000)
	w_bed = pbt.bedtool.BedTool.window_maker(pbt.BedTool(genome_file),g=genome_file,w=window_size,s=stepsize,i='srcwinnum')
	tables=[]
	tables2=[]
	for i in range(1,f+1):
		tables.append('tmp'+str(i))
		mC_bed = allc2bed('tmp'+str(i))
		mapping = pbt.bedtool.BedTool.intersect(mC_bed,w_bed,wa=True,wb=True).saveas('map'+str(i)+'.tmp')
		tables2.append('map'+str(i)+'.tmp')
		del(mapping)
	df_from_each_tmp_file = (pd.read_table(i,header=None,usecols=[10,13,6,7,8]) for i in tables2)
	m = pd.concat(df_from_each_tmp_file, ignore_index=True)
	m = m.sort_values(by = 13,ascending=True)
	f = split_df_on_column(m,size=50000000,column=13)
	c = pd.DataFrame(columns=['window','mCG_reads','CG_reads','mCG','mCHG_reads','CHG_reads','mCHG','mCHH_reads','CHH_reads','mCHH'])
	for i in range(1,f+1):
		tables.append('tmp'+str(i))
		m = pd.read_table('tmp'+str(i),header=0)
		b = window_methylation_levels(m,cutoff=cutoff,nuc_bed=(),output_mC_counts=True)
		c = pd.concat([c, b], ignore_index=True)
    if output:
        c.to_csv(output, sep='\t', index=False)
    else:
        return c
    for i in tables + tables2:
        os.remove(i)
    for i in ['w_bed','mC_bed','allc_mapping','m','b','c','df_from_each_tmp_file','tables','tables2']:
        del(i)

# count the subcontexts in fasta
def count_subcontext_fasta(fasta,context=['CG','CHG','CHH'],output=(),filter_chr=[]):
    df = pd.DataFrame(columns=['context','total_bases'])
    for c in context:
        count = 0
        for i in expand_nucleotide_code([c]):
            for sequence in SeqIO.parse(fasta, "fasta"):
                if sequence.name not in filter_chr:
                    count = count + sequence.seq.count(i) + sequence.seq.reverse_complement().count(i)
        df = df.append({'context': c, 'total_bases': count}, ignore_index=True)
    if output:
        df.to_csv(output, sep='\t', index=False)
    else:
        return df

#Count up subcontext methylation
def count_subcontext_allc(allc,context=['CG','CHG','CHH'],output=(),cutoff=3,filter_chr=[]):
    a = pd.read_table(allc)
    df = pd.DataFrame(columns=['context','total_passed','methylated','weighted_mC'])
    for c in context:
        tp = mC = t = m = 0
        for i in expand_nucleotide_code([c]):
            i_table = a[(a['mc_class']==i) & (~a['chr'].isin(filter_chr))]
            t = t + i_table['total'].values.sum()
            m = m + i_table['mc_count'].values.sum()
            i_table = i_table[i_table['total']>=cutoff]
            tp = tp + len(i_table.index)
            mC = mC + len(i_table[i_table['methylated']==1].index)
        df = df.append({'context': c, 'total_passed': tp, 'methylated': mC,
                        "weighted_mC": np.float64(m)/np.float64(t)}, ignore_index=True)
    if output:
        df.to_csv(output, sep='\t', index=False)
    else:
        return df

#Get subcontexts of genome and methylation
def subcontext_methylation(allc,fasta,context=['CG','CHG','CHH'],output=(),cutoff=3,filter_chr=[]):
    a = count_subcontext_fasta(fasta,context,output=(),filter_chr=filter_chr)
    b = count_subcontext_allc(allc,context,output=(),cutoff=3,filter_chr=filter_chr)
    df = pd.merge(a,b,on='context')
    if output:
        df.to_csv(output, sep='\t', index=False)
    else:
        return df

#output per-site methylation levels for mCs in each specified context
def per_site_mC(allc,output_path,context=['CG','CHG','CHH']):
    for i in context:
        a = filter_context(allc,[i])
        a = a[a['methylated'] == 1]
        a['mc_level'] = a['mc_count']/a['total']
        a.to_csv(output_path+i+'_site_methylation.txt', sep='\t', index=False)

#get total weighted methylation
def weighted_mC(allc, output=(), cutoff=0, genome=()):
	CG = mCG = CHG = mCHG = CHH = mCHH = CNN = mCNN = 0
	if genome:
		g = pd.read_table(genome,header=None,usecols=[0],dtype="str")
		a = filter_context(allc,['C'])
		a = a[a.chr.isin(list(g[0]))]
	else:
		a = filter_context(allc,['C'])
	for c in a.itertuples():
		if int(c[6]) >= int(cutoff):
			if c[4].startswith("CG"):
				CG = CG + int(c[6])
				mCG = mCG + int(c[5])
			elif c[4].endswith("G"):
				CHG = CHG + int(c[6])
				mCHG = mCHG + int(c[5])
			elif c[4].startswith("CN") or c[4].endswith("N"):
				CNN = CNN + int(c[6])
				mCNN = mCNN + int(c[5])
			else:
				CHH = CHH + int(c[6])
				mCHH = mCHH + int(c[5])
	b = pd.DataFrame(columns=['Context','Total','Methylated','Weighted_mC'])
	b = b.append({'Context': 'mCG', 'Total': CG, 'Methylated': mCG, 'Weighted_mC': (np.float64(mCG)/np.float64(CG))}, ignore_index=True)
	b = b.append({'Context': 'mCHG', 'Total': CHG, 'Methylated': mCHG, 'Weighted_mC': (np.float64(mCHG)/np.float64(CHG))}, ignore_index=True)
	b = b.append({'Context': 'mCHH', 'Total': CHH, 'Methylated': mCHH, 'Weighted_mC': (np.float64(mCHH)/np.float64(CHH))}, ignore_index=True)
	b = b.append({'Context': 'mCNN', 'Total': CNN, 'Methylated': mCNN, 'Weighted_mC': (np.float64(mCNN)/np.float64(CNN))}, ignore_index=True)
	if output:
		b.to_csv(output, sep='\t', index=False)
	else:
		return b
