import pandas as pd
pd.options.mode.chained_assignment = None # SettingsWithCopyWarning bug disabled , default="warn"
import argparse
import re
import functions
import helpers 
import gcsfs

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import numpy as np


def bmt_simulation_discordances(donor_vcf, host_vcf, host_sex, donor_sex, tissue, gvhd): 

	donor_annotations, numbers_dict = functions.process_annotations(donor_vcf, 'donor')
	host_annotations, host_dict = functions.process_annotations(host_vcf, 'host')

	numbers_dict.update(host_dict)
	print(donor_annotations.shape)
	print(host_annotations.shape)

	# storing the count for log entry
	donor_annotations_process_annotations = donor_annotations
	host_annotations_process_annotations = host_annotations

	func_reference = pd.read_csv("gs://mhags-data/funcotator_dataSources.v1.6.20190124g/gencode/hg19/gencode.v19.pc_transcripts_forpipeline_translated.txt", sep='\t')
	
	print('\n-------- HOST AND DONOR AFTER ANNOTATION ---------\n')
	print("\nMerging reference into donor variants")
	donor_annotations, donor_not_mutated = functions.add_funcotator_reference(donor_annotations, func_reference, 'donor')
	print("\nMerging reference into host variants")
	host_annotations, host_not_mutated = functions.add_funcotator_reference(host_annotations, func_reference, 'host')
	host_annotations = host_annotations[~pd.isna(host_annotations.mrna)]
	donor_annotations = donor_annotations[~pd.isna(donor_annotations.mrna)]

	print('\n ---------- MAKING CUSTOM PROTEOME FOR DONOR -----------')
	donor_annotations, donor_mutprots, donor_indices = functions.make_custom_proteome(donor_annotations,'donor')	
	# donor_annotations.to_csv('donor_annotations_withprot.txt' ,sep='\t', index=False)
	# donor_mutprots.to_csv('donor_mutated_proteins.txt', sep='\t', index=False)
	# donor_indices.to_csv('donor_mutation_indices.txt', sep='\t', index=False)

	print('\n ---------- MAKING CUSTOM PROTEOME FOR HOST -----------')
	host_annotations, host_mutprots, host_indices = functions.make_custom_proteome(host_annotations, 'host')
	# host_annotations.to_csv('host_annotations_withprot.txt' ,sep='\t', index=False)
	# host_mutprots.to_csv('host_mutated_proteins.txt', sep='\t', index=False)
	# host_indices.to_csv('host_mutation_indices.txt', sep='\t', index=False)

	print('\n --------- MERGING HOST AND DONOR DATA --------- \n')
	# print('Number of rows in host before merging: %s' %(host_annotations.shape[0]))
	# print("Number of genes in host annotations after removing duplicates (duplicates present due to different SNPs) %i" %(len(host_annotations["hugoSymbol"].drop_duplicates())))
	# print('Number of rows in donor before merging: %s \n' %(donor_annotations.shape[0]))

	merged = functions.merge_host_donor(host_annotations=host_annotations,
										donor_annotations=donor_annotations,
										host_mutprots=host_mutprots,
										donor_mutprots=donor_mutprots,
										host_indices=host_indices,
										donor_indices=donor_indices)

	# print('Number of rows after merging donor and host mutations: %i \n' %merged.shape[0])
	# print('\n AFTER MERGING')
	# print("\nNumber of rows after merging donor and host mutations: %i" %(merged.shape[0]))
	print(merged.shape)

	mutcount_cols = ['hugoSymbol','CHROM','start','end','refAllele','tumorSeqAllele1','tumorSeqAllele2']
	nmuts = merged[mutcount_cols].drop_duplicates()
	# print("Number of rows after merging donor and host mutations and removing duplicates: %i \n" %(nmuts.shape[0]))
	numbers_dict['Merged_mutations'] = nmuts.shape[0]

	# WRITE_PROTEOME
	print('\n --------- WRITING DONOR AND HOST PROTEOMES --------- \n')
	donor_annotations_write = merged[merged.donor_GT != '0/0']
	print('\nRows in donor after removing 0/0 genotype: %i' %donor_annotations_write.shape[0])
	donor_proteome_filename = functions.write_custom_proteome(donor_annotations_write, donor_not_mutated, 'donor')

	host_annotations_write = merged[merged.host_GT != '0/0']
	print('Rows in host after removing 0/0 genotype: %i' %host_annotations_write.shape[0])
	host_proteome_filename = functions.write_custom_proteome(host_annotations_write, host_not_mutated, 'host')

	# DO_BMT 
	print('\n -------- DOING BMT SIMULATION --------\n')
	print('\n Number of discordances before bmt simulation (should be same as merge+no duplicates): %i' %merged.shape[0])
	discordances = functions.do_bmt(discordances=merged)
	print('\n -------- DOING BMT SIMULATION --------')
	print('Number of rows in dataframe after bmt simulation: %i \n' %discordances.shape[0])

	# discordances.to_csv('discordances_after_bmt_simulation.txt', sep='\t', index=False)
	print("\nNumber of rows discordances after bmt simulation: %i" %(discordances.shape[0]))
	nmuts = discordances[mutcount_cols].drop_duplicates()
	print("Number of mutations after bmt simulation: %i" %(nmuts.shape[0]))
	numbers_dict['Number of mutations after BMT simulation'] = nmuts.shape[0]

	### for all the mutations we are interested in, find the corresponding peptide
	print('\n -------- EXTRACTING PADDED MUTATION SEQUENCES ---------\n')
	discordances = functions.extract_peptides(data=discordances, seqpad=10)
	print('\n Verify that we have the same number of lines (rows) as above: %i' %discordances.shape[0])

	print('\n ~~~~~~~~~~~ checking for NAs in the key rna and dna fields ~~~~~~~~~~~')
	# check that all the peptides have been filled in properly
	discordances_dna = discordances[pd.isnull(discordances.host_seqpad2)]
	discordances_rna = discordances[pd.isnull(discordances.host_seqpad1)]

	print('\n -------- CHOOSE THE PEPTIDE THAT IS FOREIGN TO THE DONOR ---------\n')
	discordances = functions.find_correct_sequence(data=discordances)
	print('\n Verify that we have the same number of lines (rows) as above: %i' %discordances.shape[0])

	discordances = discordances.assign(transid=discordances.annotationTranscript+':'+discordances.proteinChange)
	discordances.to_csv('all_discordances.txt', sep='\t', index=False) 

	# for printing purpose
	all_discordances = discordances
	all_discordances["start"] = all_discordances["start"].astype(int)
	all_discs_variant_classification = all_discordances["variantType"].value_counts()

	# GvHD 
	if gvhd == "yes":
		discordances_gvhd, gvhd_stats_dict = functions.gvhd_gene_filter(discordances=discordances)
		#discordances_gvhd.to_csv("gvhd_discordances.txt", sep="\t", index=False)

		print('\n -------- CHOPPING PEPTIDES INTO 8-11 SEGMENTS FOR SEARCH AND PREDICTION ---------')
		extra_cols = ['hugoSymbol','annotationTranscript','chromosome','start','end','proteinChange','variantType']

		discordances_gvhd = functions.get_all_peptides(data=discordances_gvhd, 
													mut_col='host_sequence',
													ref_col='donor_sequence', 
													othercols=extra_cols)

		discordances_gvhd = discordances_gvhd.assign(peptide_fasta_id=discordances_gvhd.hugoSymbol+'|'+
															discordances_gvhd.annotationTranscript+'|'+
															discordances_gvhd.proteinChange+'|len_'+
															discordances_gvhd.pep_length+'|num_'+
															discordances_gvhd.pep_num)
															
		discordances_gvhd = discordances_gvhd[['CHROM','start', 'end', 'tumorSeqAllele1', 'hugoSymbol','gvhd_classfication','cDnaChange','host_REF','host_ALT','donor_REF','donor_ALT','refAllele',
		'tumorSeqAllele2','variantType', 'annotationTranscript','transid','host_sequence', 'donor_sequence','proteinChange', 'seq_mut_pos', 'pep','ref_pep','pep_length','peptide_fasta_id']]
		discordances_gvhd_variantType = discordances_gvhd[["hugoSymbol", "variantType"]].drop_duplicates()
		discordances_gvhd_variantType = discordances_gvhd_variantType["variantType"].value_counts()
		discordances_gvhd["start"] = discordances_gvhd["start"].astype(int)
		discordances_gvhd.to_csv('GvHD_specific_discordances.txt', sep='\t', index=False)

		final_dict_gvhd = {"GvHD_variants_tot_num":[discordances_gvhd_variantType.to_dict()], 
		"GvHD_variants_genes":len(set(discordances_gvhd.hugoSymbol)), "GvHD_peptides":len(set(discordances_gvhd.pep))}
		final_dict_df_gvhd = pd.DataFrame(final_dict_gvhd, index=[0])
		final_dict_df_gvhd.to_csv("GvHD_dict_nums.csv", sep=",", index=False)
	
	if gvhd == "no":
		print("GvHD run not selected!")

	# GvL 
	print('\n -------- SELECTING DISCORDANCES IN TISSUE OF INTEREST ---------\n')
	tpms = functions.read_rna(param_dict['tpm'])

	discordances, genes = functions.load_gene_sets(discordances=discordances,
										tpms=tpms,   # will be none in this case used 
										just_gene_set=param_dict['just_gene_set'], # will be no which is the default 
										tissue=tissue, # will provide in terminal 
										tpms_present = param_dict['tpms_present'], # will be FALSE according to default
										sex=host_sex) # will provide in terminal 

	if param_dict['tpm'] != 'none':
			discordances = discordances.merge(tpms, left_on='hugoSymbol', right_on='Description', how='left')

	discordances = discordances[discordances.hugoSymbol.isin(genes)]
	gvl_filtered_discordances = discordances # for printing purpose 
	gvl_filtered_discordances_variantType = gvl_filtered_discordances["variantType"].value_counts()

	print("-------- SELECTING DISCORDANCES IN TISSUE OF INTEREST/SAMPLE: %s ---------" %tissue)
	print("Number of rows:discordances after filtering genes: %i" %(discordances.shape[0]))

	nmuts = discordances[mutcount_cols].drop_duplicates()
	print("Number of mutations after bmt simulation and filtering genes: %i" %(nmuts.shape[0]))
	numbers_dict['discordances_tissuespec'] = nmuts.shape[0] 
	numbers_dict['discordances_genes'] = len(set(discordances.hugoSymbol)) 
	print("Testing discordances %s" %discordances.shape[0])

	print('\n -------- CHOPPING PEPTIDES INTO 8-11 SEGMENTS FOR SEARCH AND PREDICTION ---------')
	extra_cols = ['hugoSymbol','annotationTranscript','chromosome','start','end','proteinChange','variantType']
	discordances = functions.get_all_peptides(data=discordances, 
													mut_col='host_sequence',
													ref_col='donor_sequence', 
													othercols=extra_cols)
	
	print("Total number of peptides: %i" %(len(set(discordances.pep))))
	numbers_dict['total_num_peptides'] = len(set(discordances.pep))

	## make an id for the peptides after chopping, to use for the search
	discordances = discordances.assign(peptide_fasta_id=discordances.hugoSymbol+'|'+
															discordances.annotationTranscript+'|'+
															discordances.proteinChange+'|len_'+
															discordances.pep_length+'|num_'+
															discordances.pep_num)

	discordances = discordances[['CHROM','start', 'end', 'tumorSeqAllele1','hugoSymbol','cDnaChange','host_REF','host_ALT','donor_REF','donor_ALT','refAllele','tumorSeqAllele2',
    							'variantType', 'annotationTranscript', 'host_sequence', 'donor_sequence', 'transid','proteinChange', 'seq_mut_pos', 'pep','ref_pep','pep_length','peptide_fasta_id', tissue, 'BLOOD', 'patient_exp']]

	
	discordances.to_csv(param_dict['tissue']+'_specific_discordances.txt', sep='\t', index=False)
	print(numbers_dict) # CAN BE REMOVED - JUST FOR TESTING
	peps_discordances = discordances # just for priting purpose 

	final_dict_gvl = {"Host_sex":host_sex, "Donor_sex":donor_sex, 
	"host unique genes":len(set(host_annotations["hugoSymbol"])), "donor unique genes":len(set(donor_annotations["hugoSymbol"])),
	"host_non_syn_variants":host_annotations_process_annotations.shape[0], "donor_non_syn_variants":donor_annotations_process_annotations.shape[0],
	"discordant_variants_tot_num":[all_discs_variant_classification.to_dict()], "discordant_variants_genes":len(set(all_discordances["hugoSymbol"])), 
	"GvL_variants_tot_num":[gvl_filtered_discordances_variantType.to_dict()], "GvL_variants_genes":len(set(gvl_filtered_discordances["hugoSymbol"])), "GvL_peptides":len(set(peps_discordances.pep))}
	final_dict_df_gvl = pd.DataFrame(final_dict_gvl, index=[0])
	final_dict_df_gvl.to_csv("GvL_dict_nums.csv", sep=",", index=False)

	# Known minors check
	known_minors_file = pd.read_csv("gs://mhags-data/known_minors_hg19coords_peptide.txt", sep = "\t")

	# Donor Funcotator	
	temp = pd.merge(donor_annotations_process_annotations, known_minors_file, left_on=["POS", "REF", "ALT"], right_on=["hg19_coord", "ref", "alt"])
	known_minors_file["donor_funcotator"] = known_minors_file.hg19_coord.apply(lambda x: int(x in list(temp["POS"])))

	# Host Funcotator
	temp = pd.merge(host_annotations_process_annotations, known_minors_file, left_on=["POS", "REF", "ALT"], right_on=["hg19_coord", "ref", "alt"])
	known_minors_file["host_funcotator"] = known_minors_file.hg19_coord.apply(lambda x: int(x in list(temp["POS"])))

	# BMT simulation
	temp = pd.merge(all_discordances, known_minors_file, left_on=["start", "host_REF", "host_ALT"], right_on=["hg19_coord", "ref", "alt"])
	known_minors_file["bmt_known_minors"] = known_minors_file.hg19_coord.apply(lambda x: int(x in list(temp["start"])))

	# GvHD Tissue Expression
	if gvhd == "yes":
		temp = pd.merge(discordances_gvhd, known_minors_file, left_on=["start", "host_REF", "host_ALT"], right_on=["hg19_coord", "ref", "alt"])
		known_minors_file["gvhd_known_minors"] = known_minors_file.hg19_coord.apply(lambda x: int(x in list(temp["start"])))
	
		known_minors_file.to_csv("known_minors_info.csv", index = False)

	if gvhd == "no":
		known_minors_file.to_csv("known_minors_info.csv", index = False)

def main():
	func_annotations = bmt_simulation_discordances(donor_vcf = donor_vcf, host_vcf = host_vcf, host_sex = host_sex, donor_sex = donor_sex, tissue = tissue, gvhd = gvhd) ### to do::: update the changed function name here 

if __name__ == '__main__':

	param_dict = {}

	parser = argparse.ArgumentParser()
	parser._action_groups.pop()

	parser.add_argument("-donor_vcf", help='Donor deepvariant VCF file',default=None)
	parser.add_argument("-host_vcf", help="Host deepvariant VCF file", default=None)

	parser.add_argument('-tpm', default='none', help='file path to expression tpms for tumor. name of expression column must be TPM', type=str)
	parser.add_argument('-tissue', default='AML', help='either AML, BLOOD, CML, CLL, SKIN, LUNG, COLON, default AML', type=str)
	parser.add_argument('-host_sex', default='Female', help='Male/male/M/m/ or Female/female/F/f', type=str)
	parser.add_argument('-donor_sex', default='Female', help='Male/male/M/m or Female/female/F/f', type=str)
	parser.add_argument('-use_patient_genes', default='yes', help='whether or not to use genes highly expressed in patient in the output', type=str)
	parser.add_argument('-just_gene_set', default='no', help='whether or not to just consider the gene set specified', type=str)

	#parser.add_argument("-summary_file")
	parser.add_argument("-gvhd", default = "yes", help = "To apply gvhd filters or not: only apply for AML, CML, CLL samples", type=str)

	args = parser.parse_args()
	param_dict['host_vcf'] = args.host_vcf
	param_dict['donor_vcf'] = args.donor_vcf

	param_dict['tissue'] = args.tissue
	param_dict['tpm'] = args.tpm
	param_dict['just_gene_set'] = args.just_gene_set

	param_dict['tpms_present'] = False

	if param_dict['tpm'] != 'none': param_dict['tpms_present'] = True

	if param_dict['tpms_present'] == False: param_dict['use_patient_genes'] = False
	elif args.use_patient_genes == 'yes': param_dict['use_patient_genes'] = True
	elif args.use_patient_genes == 'no': param_dict['use_patient_genes'] = False
	else: assert(False)
	
	if args.host_sex in ['Male','male','m','M']: param_dict['host_sex'] = 'M'
	elif args.host_sex in ['Female','female','f','F']: param_dict['host_sex'] = 'F'
	else: print('invalid host sex'); assert(False)

	if args.donor_sex in ['Male','male','m','M']: param_dict['donor_sex'] = 'M'
	elif args.donor_sex in ['Female','female','f','F']: param_dict['donor_sex'] = 'F'
	else: print('invalid donor sex'); assert(False)

	if (param_dict['host_sex']=='M') & (param_dict['donor_sex']=='F'): param_dict['female_to_male'] = True
	else: param_dict['female_to_male'] = False

	donor_vcf = args.donor_vcf
	host_vcf = args.host_vcf
	donor_sex = args.donor_sex
	host_sex = args.host_sex
	tissue = args.tissue
	#summary_file=args.summary_file
	gvhd = args.gvhd

	main()
	


