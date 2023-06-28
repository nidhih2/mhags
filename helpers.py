import pandas as pd
pd.options.mode.chained_assignment = None
import argparse
#import io
#import os
from Bio.bgzf import BgzfReader
#from subprocess import run, PIPE, STDOUT, Popen, check_output
import re
#import numpy as np
#import firecloud.api as fapi
#import json
#import time

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


slc=["I","L","V","F","M","C","A","G","P","T","S","Y","W","Q","N","H","E","D","K","R","Stop"]

codon=[["ATT", "ATC", "ATA"],["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"],["GTT", "GTC", "GTA", "GTG"],["TTT", "TTC"],["ATG"],["TGT", "TGC"],
		["GCT", "GCC", "GCA", "GCG"],["GGT", "GGC", "GGA", "GGG"],["CCT", "CCC", "CCA", "CCG"],["ACT", "ACC", "ACA", "ACG"],
		["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],["TAT", "TAC"],["TGG"],["CAA", "CAG"],["AAT", "AAC"],["CAT", "CAC"],["GAA", "GAG"],
		["GAT", "GAC"],["AAA", "AAG"],["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],["TAA", "TAG", "TGA"]]


aa_dict = {}
for i in range(len(codon)):
	for cd in codon[i]:
		aa_dict[cd] = slc[i]



def read_vcf_single(filename, sample):
	print("begin loading .vcf file...")
	f = BgzfReader(filename, 'r')
	lines = []
	vcf_headers = []
	colnames = []

	for line in f:
		lines.append(line[:-1].split('\t'))

	print('done reading all lines')

	for line in lines:
		if line[0].startswith('##'):
			vcf_headers.append(line[0])
		else: break


	print('found line belonging to the header')

	for i in range(len(lines)):
		if lines[i][0].startswith('#CHROM'):
			colnames = lines[i]
			body = lines[i+1:]
			break

	colnames = ['CHROM']+colnames[1:-1]+[sample]

	lines_df = pd.DataFrame(body, columns = colnames)

	print("vcf loading complete")

	return lines_df, vcf_headers

def translate(sequence):
	start_codon = 0
	end_index = len(sequence)

	if sequence.find('ATG') != -1:
		start_codon = sequence.index('ATG')
	elif sequence.find('GTG') != -1:
		start_codon = sequence.index('GTG')
		# print('\nno ATG was found, but found a GTG')
		# print(sequence)
	else: 
		print('\nno start codon was found')
		return('', '')
	# print('index of start codon: %i' %start_codon)
	proteinseq = 'M'
	for i in range(start_codon+3, len(sequence), 3):
		if len(sequence[i:i+3]) < 3: 
			proteinseq += '*'
			# print('no stop codon found')
			return(proteinseq, sequence[start_codon:])
		if aa_dict[sequence[i:i+3]] == 'Stop':
			end_index = i+3
			proteinseq += '*'
			return(proteinseq, sequence[start_codon:end_index])
		proteinseq += aa_dict[sequence[i:i+3]]

	mrna = sequence[start_codon:end_index]
	return(proteinseq, mrna)

def bmt_simulation(data, host_gt, donor_gt):

	data.reset_index(inplace=True)
	data = data.drop(labels='index', axis=1)

	datatemp = data.copy()
	hostsplit = datatemp[host_gt].str.split('/|\\|', expand=True)
	# print(hostsplit)
	donorsplit = datatemp[donor_gt].str.split('/|\\|', expand=True)
	# print(donorsplit)
	datatemp['Host_1'] = hostsplit[0]
	datatemp['Host_2'] = hostsplit[1]
	datatemp['Donor_1'] = donorsplit[0]
	datatemp['Donor_2'] = donorsplit[1]

	datatemp = datatemp[((datatemp['Host_1']!=datatemp['Donor_1']) & (datatemp['Host_1']!=datatemp['Donor_2'])) | 
				((datatemp['Host_2']!=datatemp['Donor_1']) & (datatemp['Host_2']!=datatemp['Donor_2']))]

	data = data.iloc[datatemp.index,:]
	data = data.reset_index().drop(columns='index')

	return(data)


def map_seqpad(row, sample):
	if (row.host_whichdf != 2) & (row.ALLELE=='1'): return(sample+'_seqpad1')
	if (row.host_whichdf == 2) & (row.ALLELE=='1'): return(sample+'_seqpad2')
	if (row.host_whichdf != 2) & (row.ALLELE=='0'): return(sample+'_seqpad2')
	if (row.host_whichdf == 2) & (row.ALLELE=='0'): return(sample+'_seqpad1')

	print('Not of above conditions satisfies, should be false')
	#print(row[['host_whichdf','ALLELE']])
	assert(False)
	#return(None)

def map_categ(row):
	if (row.host_whichdf != 2) & (row.ALLELE=='1'): return('c1')
	if (row.host_whichdf == 2) & (row.ALLELE=='1'): return('c2')
	if (row.host_whichdf != 2) & (row.ALLELE=='0'): return('c3')
	if (row.host_whichdf == 2) & (row.ALLELE=='0'): return('c4')
	print('not of above conditions satisfies, should be false')
	assert(False)
	return(None)

def get_gene_set(tissue, sex):

	"""
	Function for reading in gene filters, depending on the Tissue parameter provided

	@param tissue: the tissue of choice,currently options are AML, CML, CLL, BLOOD
	@param sex: either male or female, has impact on gene sets included

	@return: a list of genes

	"""
	print(sex)
	print("+++++++++++++++++++++")
	if tissue == 'AML':

		if sex in ['M', "Male"]:
			beat_male = pd.read_csv('gs://mhags-data/tissue-specific-gene-sets/FINAL_beat_aml_specific_genes_male_20210623.txt', sep='\t') # file 1 - GvL_final/AML/Filter gene set/
			genes = list(set(beat_male['symbol']))
			
		elif sex in ['F', "Female"]:
			# beat_female = pd.read_csv('gs://mhags-data/tissue-specific-gene-sets/FINAL_beat_aml_specific_genes_female_20210623.txt', sep='\t') # file 2 - GvL_final/AML/Filter gene set/
			print("This patient is female")
			beat_female = pd.read_csv('~/Desktop/beat_aml_nidhi_copy.txt', sep='\t')
			genes = list(set(beat_female['symbol']))

		elif sex == 'unknown':			
			beat_female = pd.read_csv('gs://mhags-data/tissue-specific-gene-sets/FINAL_beat_aml_specific_genes_female_20210623.txt', sep='\t') # GvL_final/AML/Filter gene set/ - REPEAT AS FILE 2
			beat_male = pd.read_csv('gs://mhags-data/tissue-specific-gene-sets/FINAL_beat_aml_specific_genes_male_20210623.txt', sep='\t') # GvL_final/AML/Filter gene set/ - REPEAT AS FILE 1
			genes = list(set(beat_female['symbol']) | set(beat_male['symbol'])) #### DOUBLE CHECK WITH KARI AS IT WAS ORIGINALLY cll_female['symbol'] instead of beat_female...
		else:
			print('Check the sex entry of donor/recipient : Invalid sex entry')

	elif tissue == "CLL":
		if sex in ['M', "Male"]:
			cll_male = pd.read_csv("gs://mhags-data/tissue-specific-gene-sets/FINAL_cll_bulk_sc_specific_genes_male_20210629.txt", sep="\t")
			genes = list(cll_male["symbol"])
		elif sex in ['F', "Female"]:
			cll_female = pd.read_csv("gs://mhags-data/tissue-specific-gene-sets/FINAL_cll_bulk_sc_specific_genes_female_20210629.txt", sep="\t")
			genes = list(cll_female["symbol"])
		elif sex == "unknown":
			cll_male = pd.read_csv("gs://mhags-data/tissue-specific-gene-sets/FINAL_cll_bulk_sc_specific_genes_male_20210629.txt", sep="\t")
			cll_female = pd.read_csv("gs://mhags-data/tissue-specific-gene-sets/FINAL_cll_bulk_sc_specific_genes_female_20210629.txt", sep="\t")
			genes = list(set(cll_male["symbol"]) | set(cll_female["symbol"]))
		else:
			print("Check the sex entry of donor/recipient : Invalid sex entry")

	elif tissue == "CML":
		if sex in ['M', "Male"]:
			cml_male = pd.read_csv("gs://mhags-data/tissue-specific-gene-sets/male_bulk_sc_combined_genes_20210623.txt", sep="\t")
			genes = list(cml_male["symbol"])
		elif sex in ['F', "Female"]:
			cml_female = pd.read_csv("gs://mhags-data/tissue-specific-gene-sets/female_bulk_sc_combined_genes_20210623.txt", sep="\t")
			genes = list(cml_female["symbol"])
		elif sex == "unknown":
			cml_m_f = pd.read_csv("gs://mhags-data/tissue-specific-gene-sets/female_bulk_sc_combined_genes_20210623.txt", sep="\t")
			genes = list(cml_m_f["symbol"])
		else:
			print("Check the sex entry of donor/recipient : Invalid sex entry")
		
	return(genes)


#### no need to reload files into the workspace as the files relating to BLOOD are already present 
def load_blood_genes(sex):

	"""
	Function to lead the blood genes

	@param sex: male or female

	@return: list of blood-specific genes

	"""

	if sex in ['M', "Male"]:
		hemato_male = pd.read_csv('gs://mhags-data/tissue-specific-gene-sets/FINAL_hemato_male_filter.txt', sep='\t') # file 1 - GvL_final/Hematopoietic Filter/ - SAME AS FILE 9 
		genes = list(hemato_male['gene_short_name'])
	elif sex in ['F', "Female"]:
		hemato_female = pd.read_csv('gs://mhags-data/tissue-specific-gene-sets/FINAL_hemato_female_filter.txt', sep='\t') # file 2 - GvL_final/Hematopoietic Filter/ - SAME AS FILE 10
		genes = list(hemato_female['gene_short_name'])
	elif sex == 'unknown':			
		hemato_m_f = pd.read_csv('gs://mhags-data/tissue-specific-gene-sets/FINAL_hemato_mf_filter.txt', sep='\t') # file 3 - GvL_final/Hematopoietic Filter/ - SAME AS FILE 11 
		genes = list(set(hemato_m_f['gene_short_name']))
	else:
		print('invalid sex entry')

	return(genes)

# 4 files
def select_patient_genes(tpms, host_sex):

	## All the tissues wrt GTEX 

	print("\n Processing Patient Gene Expression")
	print(tpms)
	tpm_above_2 = tpms[tpms['TPM'] > 2]['Description']
	print('number of tumor genes expressed above 2 TPM: %f' %len(tpm_above_2))
	print("\nNumber of genes expressed above 2 TPM: %s" %(len(tpm_above_2)))
	tpm_above_32 = tpms[tpms['TPM'] > 32]['Description']
	print('echo \"Number of genes expressed above 32 TPM: %s\" >> %s' %(len(tpm_above_32),))
 
	gtex_fem_low = pd.read_csv('gtex_rna_female_max8med5.txt', sep='\t') # FILE 1 - ~/Dropbox (Partners HealthCare)/mHAgs/mHAg_final_folder/pipeline_reference_files/gtex_summary_files/ 
	gtex_male_low = pd.read_csv('gtex_rna_male_max8med5.txt', sep='\t') # FILE 2 - gtex_summary_files/

	gprot_fem_low = pd.read_csv('gtex_prot_female_ts2.5.txt', sep='\t') # FILE 3 - gtex_summary_files/ 
	gprot_male_low = pd.read_csv('gtex_prot_male_ts2.5.txt', sep='\t') # FILE 4 - gtex_summary_files/
	
	# ts_all = pd.read_csv('')
	if host_sex in ['M', "Male"]:
		ts_low_genes = set(gprot_male_low['hgnc_symbol'])
		rna_low_genes = set(gtex_male_low['Description'])
	elif host_sex in ['F', "Female"]:
		ts_low_genes = set(gprot_fem_low['hgnc_symbol'])
		rna_low_genes = set(gtex_fem_low['Description'])
	else:
		ts_low_genes = set(gprot_fem_low['hgnc_symbol']) | set(gprot_male_low['hgnc_symbol'])
		rna_low_genes = set(gtex_male_low['Description']) | set(gtex_fem_low['Description'])

	print('nubmer TS low genes: %i' %len(ts_low_genes))

	ts_validated = set(tpm_above_2) & ts_low_genes & rna_low_genes
	print('number TS validated: %i' %len(ts_validated))
	print('echo \"Genes above 2 in patient and lowly expressed in RNA and protein GTEx data: %s\" >> %s' %(len(ts_validated)))

	other_patspec = set(tpm_above_32) & rna_low_genes
	print('number above 32 validated: %i' %len(other_patspec))
	print('echo \"Genes above 32 and lowly expressed in RNA GTEx data: %s\" >> %s' %(len(other_patspec)))

	patient_genes = list(ts_validated | other_patspec)
	print('number of patient-specific genes: %i' %len(patient_genes))
	print('echo \"Total number of patient-speciic genes: %s\" >> %s' %(len(patient_genes)))

	return(patient_genes)


def chop_peptides_basic(data, mut_col, ref_col, other_cols):

	MAX_PEP_LEN = 11
	MIN_PEP_LEN = 8

	peptides = []

	print('chopping peptides: columns')
	print(data.columns)
	

	for index, row in data.iterrows():
		mut_seq = row.at[mut_col]
		ref_seq = row.at[ref_col]

		## get length of mutant and reference sequences
		mutlen = len(mut_seq)
		reflen = len(ref_seq)

		## for most SNPs, insertions, and mutations, the mutions should be at index 10
		seq_mut_pos = 10
		## if the length of the peptide is less than 21, then it's because the mutation was close to the start or end of the protein

		## if the length of the peptide is more than 21, it's because there was an insertion, or a frameshift mutation (but the mutation should still be at location 10)
		if mutlen < 21:
			print('seqpad length less than 21')
			print('host: %s' %mut_seq)
			print('donor: %s' %ref_seq)
			for k in range(min(mutlen, reflen)):
				if mut_seq[k] != ref_seq[k]: 
					seq_mut_pos = k
					break

		for pep_length in range(MIN_PEP_LEN, MAX_PEP_LEN+1):
			# print('peptide length: %i' %pep_length)

			first_start_index = int(max(0, seq_mut_pos-pep_length+1))
			last_start_index = int(max(mutlen, reflen))-pep_length+1

			for i in range(first_start_index, last_start_index):
				mut_smolpep = mut_seq[i:i+pep_length]
				ref_smolpep = ref_seq[i:i+pep_length]

				pep_mutloc = seq_mut_pos - (pep_length)

				peptides.append([row.at[mut_col], row.at[ref_col], mut_smolpep, ref_smolpep, str(i), pep_mutloc, str(pep_length), seq_mut_pos+1]+list(data.loc[index, other_cols]))

				

	pepdf = pd.DataFrame(peptides, columns=['host_sequence', 'donor_sequence', 'pep', 'ref_pep', 'pep_num', 'pep_mut_loc', 'pep_length', 'seq_mut_pos']+other_cols)

	# pepdf = pepdf.drop(columns='index')
	pepdf = pepdf.drop_duplicates()

	just_mut = set(pepdf.pep) - set(pepdf.ref_pep)
	print('number of peptides unique to host: %i' %len(just_mut))

	pepdf = pepdf[~pepdf.pep.isin(pepdf.ref_pep)]
	# pepdf = pepdf[pepdf.pep.isin(just_mut)]
	print('number of peptides for hlathena: %i' %pepdf.shape[0])

	return(pepdf)

