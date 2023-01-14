import pandas as pd
pd.options.mode.chained_assignment = None
import re
import helpers as helpers
import translation_helpers as thelp
import concurrent.futures
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import openpyxl

colsp = ['hugoSymbol','chromosome','transcriptStrand','variantClassification','variantType','refAllele','tumorSeqAllele1','tumorSeqAllele2','annotationTranscript',
		'cDnaChange','codonChange','proteinChange']

cols_keep = ['hugoSymbol', 'ncbiBuild', 
			'chromosome', 'start', 'end', 
			'variantClassification', 'secondaryVariantClassification', 'variantType', 
			'refAllele', 'tumorSeqAllele1', 'tumorSeqAllele2', 
			'genomeChange', 'annotationTranscript', 'transcriptStrand', 'transcriptExon',
			'transcriptPos', 'cDnaChange', 'codonChange', 'proteinChange', 
			'gcContent', 'referenceContext', 'otherTranscripts']


def process_annotations(vcf_filename, sample):
	# annotations = annotations[annotations.FILTER=='PASS'] ## can remove when filtering is properly in place

	temp_numbers_logfile = {}

	annotations, header = helpers.read_vcf_single(vcf_filename, sample)
	print(annotations.shape)

	ann_colnames = ''
	for head in header:
		if head.startswith('##INFO=<ID=FUNCOTATION'):
			ann_colnames=head
	# print(ann_colnames)

	print('Number of rows in %s funcotations: %i' %(sample, annotations.shape[0])) # 24868

	#temp_numbers_logfile[sample+'_vcf_rows'] = annotations.shape[0]

	#print(ann_colnames)
	ann_colnames = ann_colnames.replace('\"', '').replace('>','').split(': ')
	# print(ann_colnames)
	ann_colnames = ann_colnames[1].split('|')
	ann_colnames = pd.Series(ann_colnames).str.replace('Gencode_19_','')

	annotations.INFO = annotations.INFO.str.replace('FUNCOTATION=', '')
	annotations.INFO = annotations.INFO.str.split('#')
	#annotations.to_csv("annotations_to_compare.txt", index=False, sep="\t")
	annotations = annotations.explode('INFO')
	#annotations.to_csv("donor_annotations_explode.txt", index=False, sep="\t")
	annotations.INFO = annotations.INFO.str.split(',')
	print(annotations.shape)
	annotations = annotations.explode('INFO')
	print(annotations.shape)
	annotations.INFO = annotations.INFO.str.replace(r'\[|\]', '', regex=True)
	annotations = annotations[~ annotations.INFO.str.startswith('Unknown', na=False)]
	muts = annotations.INFO.str.split(r'\|', expand=True)
	muts.columns = ann_colnames
	muts = muts[cols_keep]

	annotations = pd.concat([annotations, muts], axis=1)
	annotations = annotations[annotations.proteinChange!='']

	## remove any annotations that we don't care about
	annotations = annotations[~annotations.variantClassification.isin(['IGR','INTRON','FIVE_PRIME_FLANK','THREE_PRIME_UTR','FIVE_PRIME_UTR','SPLICE_SITE','COULD_NOT_DETERMINE','RNA','LINCRNA','NONSTOP','SILENT'])]
	annotations = annotations[annotations.variantType!='ONP']

	print('number of rows in %s funcotations after removing IGR and other unwanted mutation types: %i' %(sample, annotations.shape[0]))
	temp_numbers_logfile[sample+'_varClass_filter'] = annotations.shape[0]

	## split the genotype field so we can see the genotype
	annotations[sample+'_GT'] = annotations[sample].str.split(':', expand=True)[0]

	print(annotations[sample+'_GT'].value_counts())

	## remove any loci that for some reason are 0/0
	annotations = annotations[annotations[sample+'_GT'].isin(['0/1','1/1','1/2'])]

	print('\nNumber of rows in %s funcotations after removing loci other than \'0/1\',\'1/1\',\'1/2\': %i' %(sample, annotations.shape[0]))
	temp_numbers_logfile[sample+'_GT_filter'] = annotations.shape[0]

	loci_add = {"0/1":'het', "1/2":'het', "1/1":'hom'}
	annotations[sample+'_loci'] = annotations[sample+'_GT'].map(loci_add)

	print('\n Number of transcripts in annotations: %i' %len(set(annotations.annotationTranscript)))
	print('\nMutation breakdown:')
	print(annotations.variantClassification.value_counts())
	print('\nType breakdown:')
	print(annotations.variantType.value_counts())

	return(annotations, temp_numbers_logfile)


def add_funcotator_reference(annotations, func_reference, sample):

	"""
	merge in the reference sequence data from the Funcotator reference funcotator_dataSources.v1.6.20190124g/gencode/hg19

	@param annotations_filename: vcf file location as output from get_func_results
	@param func_reference: the funcotator reference sequences to merge
	@param sample: either host or donor

	@return: the merged annotations, sequences from the reference that have no mutations in the sample
	"""

	## remove rows of the func reference whose mRNA sequence is null in the row 
	#print(annotations)
	func_reference = func_reference[~pd.isnull(func_reference.mrna)]

	print('Number of rows in %s vcf file: %i' %(sample, annotations.shape[0]))
	annotations = pd.merge(annotations, func_reference, left_on='annotationTranscript', right_on='ensembl_transcript', how='inner').reset_index().drop('index', axis=1)
	print('Number of rows in %s vcf file after merging in reference: %i' %(sample, annotations.shape[0]))

	not_mutated = func_reference[~func_reference.ensembl_transcript.isin(annotations.annotationTranscript)].reset_index().drop(columns='index')
	print('Number of transcripts with a mutation %i' %(len(set(annotations.ensembl_transcript))))
	print('Number of transcripts with no mutation %i' %(len(set(not_mutated.ensembl_transcript))))

	## a weird edge case (from AML003) that I don't have the patience to deal with at the moment
	annotations = annotations[~annotations.proteinChange.str.contains('*', regex=False) & (annotations.variantClassification!='IN_FRAME_DEL')]

	annotations = annotations[~pd.isnull(annotations.cDnaChange)].reset_index().drop('index', axis=1)

	## extract information from cDNA HGVS - splitting the cDnachange into 4 columns of all the annotations and integrating with the annotations 
	parsed = thelp.parse_cdna_change(annotations, 'cDnaChange')
	annotations = pd.concat([annotations, parsed], axis=1)

	annotations = annotations[~pd.isnull(annotations.proteinChange)].reset_index().drop('index', axis=1)

	## extract information from protein HGVS
	prot_parsed = thelp.parse_protein_change(annotations)
	annotations = pd.concat([annotations, prot_parsed], axis=1)

	return annotations, not_mutated


def custom_proteome_thread_function(name, genes_to_mutate, annotations, sample):
	"""
	Function for applying mutations to the mrna sequences. this is run on multiple threads.

	@param name: the name of the thread (chromosome)
	@param genes_to_mutate: genes on the chromosome that need to be processed
	@param annotations: the data used for applying mutations
	@param sample: either host or donor

	@return: a list with the three pieces of information needed: the mutated seuqences, 
					indices of each mutation, 
					and a dataframe encoding which strand each mutations is on
	"""
	#logging.info("Thread %s: starting", name)

	genes_version = pd.DataFrame(columns=['annotationTranscript',sample+'_cdnamut1',sample+'_cdnamut2',sample+'_protmut1',sample+'_protmut2'])
	whichdf = pd.DataFrame(columns=['annotationTranscript','cDnaChange','proteinChange',sample+'_whichdf'])


	counter = 0
	genes_to_mutate = pd.unique(annotations.annotationTranscript)
	mut_indices = pd.DataFrame(columns=['annotationTranscript','cDnaChange','proteinChange',sample+'_base_index_1',sample+'_base_index_2'])

	# newmerge_keep = pd.DataFrame(columns=['annotationTranscript','proteinChange1','proteinChange2',sample+'_mut_index',sample+'_seqpad1',sample+'_seqpad2'])
	# newmerge_keep = pd.DataFrame()

	for gene in genes_to_mutate:
		## subset to just mutations in a transcript
		sub = annotations[annotations.annotationTranscript==gene]

		## for some reason, sorting doesn't always work because of a typing error. remove this later
		try:
			sub = sub.sort_values(by='mutstart').reset_index().drop('index', axis=1)
		except:
			print('sorting didn\'t work for some reason')
			print(sub)
			print(sub.mutstart)
			assert(False)

		## function to sort the rows/mutations in to two different tables based on which strand they should be on
		## table1 by default has the mutation, and table2 the wild type (if het), but there are special cases involving mutations
		##		that overlap, and framshift mutations are always put in the second table if possible
		sub, table1, table2 = thelp.determine_mutation_order(sub, sample)

		whichdf = pd.concat([whichdf,sub[['annotationTranscript','cDnaChange','proteinChange',sample+'_whichdf']]])

		## need a different columns for each protein strand
		mutind1 = thelp.find_mutation_indices(table1, sample,'1')
		mutind2 = thelp.find_mutation_indices(table2, sample,'2')
		mutinds = mutind1.merge(mutind2, how='outer')
		mut_indices = pd.concat([mut_indices, mutinds])

		## other option is to ignore and remove these from the dataset entirely
		if len(sub.loc[0,'refProtein']) < max(sub.protstart):
			print('the protein is way too short: %s (%i long), the first and last mutation are in position %i and %i' %(sub.loc[0,'refProtein'], len(sub.loc[0,'refProtein']), min(sub.protstart), max(sub.protstart)))
			# print('continuing')
			sub_version = sub[['annotationTranscript','cDnaChange','proteinChange']]
			sub_version[sample+'_cdnamut1'] = sub.mrna
			sub_version[sample+'_cdnamut2'] = sub.mrna
			sub_version[sample+'_protmut1'] = sub.refProtein
			sub_version[sample+'_protmut2'] = sub.refProtein
			# continue
		else:
			sub_version = thelp.extract_full_proteins(sub, table1, table2, gene, sample)

		genes_version = pd.concat([genes_version, sub_version])
		

	print("Thread for chromosome %s: finishing" %name)
	return([genes_version, mut_indices, whichdf])


def make_custom_proteome(annotations, sample):
	"""
	A function for applying all the mutations in the host or donor to the reference sequences to produce a custom proteome, and does so threading on chromosome

	@param annotations: the data containing mutations to apply and reference sequences
	@param sample: either host or donor
	@param logfile: file for logging milestone in the pipeline

	@return: the annotations with strand information merged, the mutated proteins, and the mutation indices

	"""
	# null = annotations[pd.isnull(annotations.protstart) | pd.isnull(annotations.protend)]
	to_drop = list(set(['QUAL','FILTER','INFO', 'FORMAT','ID']) & set(annotations.columns))
	annotations = annotations.drop(columns=to_drop)

	# na = annotations[pd.isnull(annotations.protstart)]
	# print(na[['hugoSymbol','annotationTranscript','variantType','variantClassification','proteinChange','cDnaChange','mutstart','mutend','mutref','mutalt','protstart','protend','protref','protalt']])
	# assert(False)

	annotations.protstart = annotations.protstart.astype(int)
	annotations.protend = annotations.protend.astype(int)
	annotations.mutstart = annotations.mutstart.astype(int)
	annotations.mutend = annotations.mutend.astype(int)

	annotations['mutref'] = annotations['mutref'].fillna('')
	annotations['mutalt'] = annotations['mutalt'].fillna('')
	annotations['protref'] = annotations['protref'].fillna('')
	annotations['protalt'] = annotations['protalt'].fillna('')

	## dataframe that holds the mutant protein sequences
	genes_version = pd.DataFrame(columns=['annotationTranscript',sample+'_cdnamut1',sample+'_cdnamut2',sample+'_protmut1',sample+'_protmut2'])
	whichdf = pd.DataFrame(columns=['annotationTranscript','cDnaChange','proteinChange',sample+'_whichdf'])

	counter = 0
	genes_to_mutate = pd.unique(annotations.annotationTranscript)
	mut_indices = pd.DataFrame(columns=['annotationTranscript','cDnaChange','proteinChange',sample+'_base_index_1',sample+'_base_index_2'])

	# newmerge_keep = pd.DataFrame(columns=['annotationTranscript','proteinChange1','proteinChange2',sample+'_mut_index',sample+'_seqpad1',sample+'_seqpad2'])
	newmerge_keep = pd.DataFrame()

	executor = concurrent.futures.ThreadPoolExecutor(max_workers=23)
	threadresults = list()
	for chrom in pd.unique(annotations.CHROM):
		ann_sub = annotations[annotations.CHROM==chrom]
		genes_to_mutate = pd.unique(ann_sub.annotationTranscript)
		print("%s: number of mutations in chromosome %s: %i" %(sample, chrom, ann_sub.shape[0]))
		x = executor.submit(custom_proteome_thread_function, chrom, genes_to_mutate, ann_sub, sample)
		threadresults.append(x)

	genes_versions = []
	mut_indiceses = []
	whichdfs = []
	for future in concurrent.futures.as_completed(threadresults):
		genes_versions.append(future.result()[0])
		mut_indiceses.append(future.result()[1])
		whichdfs.append(future.result()[2])

	genes_version = pd.concat(genes_versions).reset_index().drop(columns='index')
	mut_indices = pd.concat(mut_indiceses).reset_index().drop(columns='index')
	whichdf = pd.concat(whichdfs).reset_index().drop(columns='index')

	mut_indices = mut_indices.fillna(0)

	annotations = annotations.merge(whichdf)

	annotations = annotations.reset_index().drop(columns='index')

	return(annotations, genes_version, mut_indices)

def merge_host_donor(host_annotations, donor_annotations, host_mutprots, donor_mutprots, host_indices, donor_indices):
	
	"""
	Function for merging together the host and donor annotations, important for doing the bmt simulation later in the pipeline
	Need to make sure that any mutation that both donor and host have end up on the same line, so the simulation can be done properly

	@param host_annotations: dataframe with the host annotations
	@param donor_annotations: data frame withthe donor annotations
	@param host_mutprots: dataframe with a gene per line, with both mutated proteins for the host
	@param donor_mutprots: dataframe with a gene per line, with both mutated protein for the donor
	@param host_indices: dataframe with the index shift for the host
	@param donor_indices: dataframe with the index shift for the donor
	@param case_name: unique identifier for this pipeline run

	@return: dataframe with merged data

	"""
	host_annotations = host_annotations.rename(columns={'REF':'host_REF', 'ALT':'host_ALT'})
	donor_annotations = donor_annotations.rename(columns={'REF':'donor_REF', 'ALT':'donor_ALT'})

	## just label which sites are '1/2', meaning two mutations in the same place
	host_annotations = host_annotations.assign(host_hetsite='no')
	for index, row in host_annotations.iterrows():
		if row.host_GT=='1/2': host_annotations.loc[index, 'host_hetsite'] = 'yes'


	merge_1 = host_annotations.merge(donor_annotations, how='outer')
	#print(len(set(host_annotations.proteinChange) & set(donor_annotations.proteinChange)))
	#print(set(host_annotations.columns) & set(donor_annotations.columns))
	print("\nNumber of rows after merging host_annotations and donor_annotations: %i" %(merge_1.shape[0]))
	# assert(False)
	host_add = list(set(merge_1.annotationTranscript) - set(host_annotations.annotationTranscript))

	## find transcripts in the donor that I need to add in, not present in the host
	host_add = merge_1[merge_1.annotationTranscript.isin(host_add)][['annotationTranscript','mrna','refProtein']]
	host_add = host_add.assign(host_protmut1=host_add.refProtein)
	host_add = host_add.assign(host_protmut2=host_add.refProtein)
	host_add = host_add.assign(host_cdnamut1=host_add.mrna)
	host_add = host_add.assign(host_cdnamut2=host_add.mrna)
	host_add = host_add.drop(columns=['mrna','refProtein']).drop_duplicates()

	host_mutprots = pd.concat([host_mutprots, host_add]).reset_index().drop(columns='index')
	#print(host_mutprots[['annotationTranscript','cDnaChange','proteinChange','host_protmut1','host_protmut2']])
	host_na = host_mutprots[pd.isnull(host_mutprots.host_protmut1) | pd.isnull(host_mutprots.host_protmut2)]
	print('\n Number of lines in host that have NA for protmut: %i' %host_na.shape[0])
	#print(host_na[['annotationTranscript','cDnaChange','proteinChange','host_protmut1','host_protmut2']])

	## this is just for the test cases:
	if (('cDnaChange' in host_mutprots.columns) & ('proteinChange' in host_mutprots.columns)):
		host_mutprots = host_mutprots.drop(columns=['cDnaChange','proteinChange']).drop_duplicates().reset_index().drop(columns='index')

	#print(host_mutprots)
	#print(host_mutprots.columns)

	merge_1 = merge_1.merge(host_mutprots, how='left')
	print("\nnumber of rows after merging host mutated proteins:: %i" %(merge_1.shape[0]))

	donor_add = list(set(merge_1.annotationTranscript) - set(donor_annotations.annotationTranscript))

	## find trnascripts in the host that I need to add in, not present in the donor
	donor_add = merge_1[merge_1.annotationTranscript.isin(donor_add)][['annotationTranscript','mrna','refProtein']]
	donor_add = donor_add.assign(donor_protmut1=donor_add.refProtein)
	donor_add = donor_add.assign(donor_protmut2=donor_add.refProtein)
	donor_add = donor_add.assign(donor_cdnamut1=donor_add.mrna)
	donor_add = donor_add.assign(donor_cdnamut2=donor_add.mrna)
	donor_add = donor_add.drop(columns=['mrna','refProtein']).drop_duplicates()

	## add in reference sequences for any genes present in one patient but not the other
	donor_mutprots = pd.concat([donor_mutprots, donor_add]).reset_index().drop(columns='index')
	#print(donor_mutprots[['annotationTranscript','cDnaChange','proteinChange','donor_protmut1','donor_protmut2']])
	donor_na = donor_mutprots[pd.isnull(donor_mutprots.donor_protmut1) | pd.isnull(donor_mutprots.donor_protmut2)]
	print('\n number of lines in donor that have NA for protmut: %i' %donor_na.shape[0])
	#print(donor_na[['annotationTranscript','cDnaChange','proteinChange','donor_protmut1','donor_protmut2']])
	
	## this is just for the test cases:
	if (('cDnaChange' in donor_mutprots.columns) & ('proteinChange' in donor_mutprots.columns)):
		donor_mutprots = donor_mutprots.drop(columns=['cDnaChange','proteinChange']).drop_duplicates().reset_index().drop(columns='index')

	#print(donor_mutprots)
	#print(donor_mutprots.columns)

	merge_1 = merge_1.merge(donor_mutprots, how='left')
	print("number of rows after merging donor mutated proteins: %i" %(merge_1.shape[0]))

	merge_1 = merge_1.merge(host_indices, how='left')
	print("number of rows after merging host indices: %i" %(merge_1.shape[0]))

	merge_1 = merge_1.merge(donor_indices, how='left')
	print("number of rows after merging donor indices: %i" %(merge_1.shape[0]))

	## need to do some filling NAs such that the data is comparable
	merge_1.host_GT = merge_1.host_GT.fillna('0/0')
	merge_1.donor_GT = merge_1.donor_GT.fillna('0/0')

	for col in ['protref','protalt']:
		merge_1[col] = merge_1[col].fillna('')

	merge_1.mutstart = merge_1.mutstart.astype(int)
	merge_1.mutend = merge_1.mutend.astype(int)
	merge_1 = merge_1.sort_values('mutstart').reset_index().drop(columns='index')


	unique_muts = merge_1[['hugoSymbol','CHROM','refAllele','tumorSeqAllele1','tumorSeqAllele2']].drop_duplicates()
	print("\nnumber of unique mutations in merged file: %i" %(unique_muts.shape[0]))

	return(merge_1)


def write_custom_proteome(mutated_sequences, reference_sequences, sample):
	"""
	Function for creating a fasta file with the complete custom proteome for one of the patients

	@param mutated_sequences: dataframe with sequences that contains mutations
	@param reference_sequences: dataframe with any sequence that does not have a variant, with the reference protein
	@param folder_name: folder where data should be saved
	@param case_name: unique identifier for the pipeline run
	@param sample: host or donor

	@return: the filename of the custom proteome

	"""

	in_both = list(set(mutated_sequences.columns) & set(reference_sequences.columns))

	mut_merge = mutated_sequences[in_both+[sample+'_protmut1', sample+'_protmut2']].drop_duplicates().reset_index().drop(columns='index')

	
	ref_merge = reference_sequences[in_both]
	ref_merge[sample+'_protmut1'] = ref_merge.refProtein
	ref_merge[sample+'_protmut2'] = ref_merge.refProtein
	#print('\n reference sequences:')
	#print(ref_merge)
	ref_na = ref_merge[pd.isnull(ref_merge.mrna)]
	print('rows in ref_merge that have na for mrna: %i' %ref_na.shape[0])

	print('rows in proteome before merging: %i' %mut_merge.shape[0])

	full_proteome = pd.concat([mut_merge, ref_merge], axis=0)
	print('rows in proteome after merging: %i' %full_proteome.shape[0])

	na = full_proteome[pd.isnull(full_proteome[sample+'_protmut1']) | pd.isnull(full_proteome[sample+'_protmut2'])]

	#print(full_proteome)
	#print(na[['gene','ensembl_gene','ensembl_transcript','mrna','refProtein',sample+'_protmut1',sample+'_protmut2']])
	
	full_proteome = full_proteome[~pd.isnull(full_proteome.mrna)]

	# assert(False)

	prot_filename = sample+'_custom_proteome.fasta'
	
	with open(prot_filename, 'w') as f:
		for index, row in full_proteome.iterrows():
			f.write('>'+row.gene+'|'+row.ensembl_gene+'|'+row.ensembl_transcript+'|transcript_1\n')
			# f.write(row[sample+'_protmut1']+'\n')
			f.write(row[sample+'_protmut1']+'\n')
			f.write('>'+row.gene+'|'+row.ensembl_gene+'|'+row.ensembl_transcript+'|transcript_2\n')
			f.write(row[sample+'_protmut2']+'\n')
			# f.write(row[sample+'_protmut2']+'\n')
	f.close()

def do_bmt(discordances):

	"""
	Function for doing bone marrow simulation to find directional variants (allele found in host but not in donor)

	@param discordances: dataframe with merged data

	@return: a data frame with just variants with the correct directionality

	"""

	discordances = discordances[discordances.variantClassification!='SILENT']
	## also need to replace all 1/2 with 1/1, so that recognized peptides are actually removed
	
	discordances.host_GT = discordances.host_GT.replace('1/2', '1/1')
	discordances.donor_GT = discordances.donor_GT.replace('1/2', '1/1')

	discordances = helpers.bmt_simulation(discordances, 'host_GT', 'donor_GT')

	return(discordances)

def extract_peptides(data, seqpad):

	"""
	Function for pulling out the sequences surrounding each of the variants

	@param data: the discordant variants
	@param seqpad: the number of amino acids to include to the left and right of each mutation

	@return: the data with sequence padded variants added in

	"""

	mutpeps1 = thelp.fetch_seqpad(data, 'host', seqpad)
	mutpeps2 = thelp.fetch_seqpad(data, 'donor', seqpad)

	newmerge = mutpeps1.merge(mutpeps2, how='outer')
	data = data.merge(newmerge)
	return(data)

def find_correct_sequence(data):

	"""
	Function for picking the sequence that is foreign to the donor
	Because the reference allele can be the foreign one, we need to decide which allele is the foreign allele,
		and then choose the seuqence accordingly

	@param data: the dataframe with all the discordances
	@param logfile: file for logging progress

	@return: the data with the approprate sequence labeled as 'host_sequence'
	"""
	print('\n*** Finding Correct Sequence ***\n')

	data.host_whichdf = data.host_whichdf.fillna(0)

	host = data.host_GT.str.split(r"/|\|", expand=True)
	host.columns = ['host_GT_'+str(col) for col in host.columns]
	donor = data.donor_GT.str.split(r"/|\|", expand=True)
	donor.columns = ['donor_GT_'+str(col) for col in donor.columns]

	genos = pd.concat([host, donor], axis=1)

	alleles = []
	for index, row in genos.iterrows():
		if (row.host_GT_0 != row.donor_GT_0) & (row.host_GT_0 != row.donor_GT_1):
			alleles.append(row.host_GT_0)
		elif (row.host_GT_1 != row.donor_GT_0) & (row.host_GT_1 != row.donor_GT_1):
			alleles.append(row.host_GT_1)
		else: 
			print('neither allele is host-specific... you done')
			assert(False)

	data = data.assign(ALLELE=alleles)

	print("\nDistribution of antigenic allele:" %(data.ALLELE.value_counts()))
	#print("%s\" >> %s' %(data.ALLELE.value_counts(), logfile))

	data = data.assign(host_sequence_id=data.apply(lambda row: helpers.map_seqpad(row, 'host'), axis=1))
	data = data.assign(host_sequence=data.lookup(data.index, data.host_sequence_id))
	data = data.assign(host_condition=data.apply(lambda row: helpers.map_categ(row), axis=1))

	data = data.assign(donor_sequence_id=data.apply(lambda row: helpers.map_seqpad(row, 'donor'), axis=1))
	data = data.assign(donor_sequence=data.lookup(data.index, data.donor_sequence_id))
	return(data)

def read_rna(tpm_filename):

	"""
	A function for reading the gct formatted file with RNA tpm data

	@param tpm_filename: file path to location of TPMs
	@param logfile: file for logging progress

	@return: dataframe with TPMs

	"""
	if tpm_filename != 'none':
		tpms = pd.read_csv(tpm_filename)
		tpms = tpms.iloc[2:,]
		tpms = tpms.loc[:,'#1.2'].str.split('\\t', expand=True)
		tpms.columns = ['Name', 'Description', 'TPM']
		tpms['TPM'] = pd.to_numeric(tpms['TPM'])
		print("\nPatient RNAseq loading complete.")

	else: tpms = 'none'

	return(tpms)

def load_gene_sets(discordances, tpms, just_gene_set, tissue, tpms_present, sex):
	
	"""
	Function for leading the blood genes, and for selecting patient-specific genes in the TPMs are provided.

	@param 

	"""
	gene_set_genes = helpers.get_gene_set(tissue=tissue, sex=sex)

	blood_genes = []
	if tissue != 'BLOOD':
		## load the blood genes separately, so we can do everything in one run
		blood_genes = helpers.load_blood_genes(sex=sex)

	## FIND THE GENES THAT ARE UPREGULATED IN THE PATIENT
	if tpms_present:
		patient_genes = helpers.select_patient_genes(tpms=tpms, host_sex=sex)
	# else: patient_genes = list(discordances.hugoSymbol)
	else: patient_genes = []

	if just_gene_set == 'yes':
		patient_genes = []
		blood_genes= []

	discordances[tissue] = 'no'
	discordances.loc[discordances['hugoSymbol'].isin(gene_set_genes), tissue] = 'yes'

	discordances['BLOOD'] = 'no'
	discordances.loc[discordances['hugoSymbol'].isin(blood_genes), 'BLOOD'] = 'yes'

	discordances['patient_exp'] = 'no'
	discordances.loc[discordances['hugoSymbol'].isin(patient_genes), 'patient_exp'] = 'yes'

	allgenes = list(set(gene_set_genes) | set(blood_genes) | set(patient_genes))

	print("\nNumber of patient genes: %i" %(len(patient_genes)))
	print("Number of %s genes: %i" %(tissue, len(gene_set_genes)))
	print("Number of blood genes: %i" %(len(blood_genes)))
	print("\nNumber of genes to consider: %i" %(len(allgenes)))

	return(discordances, allgenes)

def get_all_peptides(data, mut_col, ref_col, othercols):

	"""
	Function to process data and call another function to chop the 'host_seqpad' sequence into 8-11-mers.

	@param data: data with discordances and sequences
	@param mut_col: column name to use for foreign seuqence from the host
	@param ref_col: column name to use for the ref seuqence from the donor

	@return the data with all the peptides split into 8-11-mers. Should be a large dataframe

	"""

	data[mut_col] = data[mut_col].str.replace(r'\*','', regex=True)
	data[ref_col] = data[ref_col].str.replace(r'\*','', regex=True)

	## do the chopping
	mut_peptides = helpers.chop_peptides_basic(data, mut_col, ref_col, othercols)

	mut_peptides = mut_peptides.drop_duplicates().reset_index().drop(columns='index')

	tomerge = list(set(mut_peptides.columns) & set(data.columns))
	#print('\ncolumns we are merging on:')
	#print(tomerge)

	data = data.merge(mut_peptides)
	# assert(False)
	return(data)

def gvhd_gene_filter(discordances):
    discordances = discordances
    gvhd_gene_classification = pd.read_csv("gs://mhags-data/gvhd-tissues-data/gvhd_gene_list.txt", sep="\t")

    discordances_gvhd = discordances[discordances["hugoSymbol"].isin(gvhd_gene_classification["gvhd_genes"])]
    discordances_gvhd = discordances_gvhd.merge(gvhd_gene_classification, left_on="hugoSymbol", right_on="gvhd_genes", how="inner")
    
    gvhd_classes = discordances_gvhd[["hugoSymbol", "gvhd_classfication"]].drop_duplicates()
    gvhd_classes_specific = gvhd_classes["gvhd_classfication"].value_counts().to_dict()
    
    gvhd_classes["gvhd_classfication_list"] = gvhd_classes["gvhd_classfication"].str.split(",")
    gvhd_tissue_bins = gvhd_classes["gvhd_classfication_list"].str.len().value_counts().to_dict()
    
    gvhd_num_dict = {"Pan expressed genes":gvhd_classes_specific["pan expressed gene"],
            "Genes expressed in non-hemato tissues":gvhd_classes_specific["lung,liver,skin,colon,lacrimal,oral"],
            "Genes expressed in lung":gvhd_classes_specific["lung"], 
            "Genes expressed in liver":gvhd_classes_specific["liver"], 
            "Genes expressed in skin":gvhd_classes_specific["skin"], 
            "Genes expressed in colon":gvhd_classes_specific["colon"],
            "Genes expressed in lacrimal":gvhd_classes_specific["lacrimal"],
            "Genes expressed in oral":gvhd_classes_specific["oral"],
            "Genes expressed in hematopoietic":gvhd_classes_specific["hematopoietic"],
            "Genes expressed by 2 gvhd tissues":gvhd_tissue_bins[2],
            "Genes expressed by 3 gvhd tissues":gvhd_tissue_bins[3],
            "Genes expressed by 4 gvhd tissues":gvhd_tissue_bins[4],
            "Genes expressed by 5 gvhd tissues":gvhd_tissue_bins[5]}

    return(discordances_gvhd, gvhd_num_dict)

#def gvhd_tissues(discordances):
# 	#discordances = pd.read_csv("all_discordances.txt", sep="\t")
# 	lung = pd.read_excel("gs://mhags-data/gvhd-tissues-data/specificGenes_7tissues.xlsx", sheet_name="lung")
# 	liver = pd.read_excel("gs://mhags-data/gvhd-tissues-data/specificGenes_7tissues.xlsx", sheet_name="liver")
# 	skin = pd.read_excel("gs://mhags-data/gvhd-tissues-data/specificGenes_7tissues.xlsx", sheet_name="skin")
# 	colon = pd.read_excel("gs://mhags-data/gvhd-tissues-data/specificGenes_7tissues.xlsx", sheet_name="colon")
# 	lacrimal = pd.read_excel("gs://mhags-data/gvhd-tissues-data/specificGenes_7tissues.xlsx", sheet_name="lacrimal")
# 	oral = pd.read_excel("gs://mhags-data/gvhd-tissues-data/specificGenes_7tissues.xlsx", sheet_name="oral")
# 	hemato = pd.read_excel("gs://mhags-data/gvhd-tissues-data/specificGenes_7tissues.xlsx", sheet_name="hemato")

# 	lung_genes= list(lung["genes"])
# 	liver_genes = list(liver["genes"])
# 	skin_genes = list(skin["genes"])
# 	colon_genes = list(colon["genes"])
# 	lacrimal_genes = list(lacrimal["genes"])
# 	oral_genes = list(oral["genes"])
# 	hemato_genes = list(hemato["genes"])

# 	print("The total number of lung specific genes are %i" %len(lung_genes))
# 	print("The total number of liver specific genes are %i" %len(liver_genes))
# 	print("The total number of skin specific genes are %i" %len(skin_genes))
# 	print("The total number of colon specific genes are %i" %len(colon_genes))
# 	print("The total number of lacrimal specific genes are %i" %len(lacrimal_genes))
# 	print("The total number of oral specific genes are %i" %len(oral_genes))
# 	print("The total number of hematopoeitic specific genes are %i" %len(hemato_genes))

# 	lung_specific_genes = discordances[discordances.hugoSymbol.isin(lung_genes)]["hugoSymbol"]
# 	liver_specific_genes = discordances[discordances.hugoSymbol.isin(liver_genes)]["hugoSymbol"]
# 	skin_specific_genes = discordances[discordances.hugoSymbol.isin(skin_genes)]["hugoSymbol"]
# 	colon_specific_genes = discordances[discordances.hugoSymbol.isin(colon_genes)]["hugoSymbol"]
# 	lacrimal_specific_genes = discordances[discordances.hugoSymbol.isin(lacrimal_genes)]["hugoSymbol"]
# 	oral_specific_genes = discordances[discordances.hugoSymbol.isin(oral_genes)]["hugoSymbol"]
# 	hemato_specific_genes = discordances[discordances.hugoSymbol.isin(hemato_genes)]["hugoSymbol"]

# 	discordances['gvhd_filter'] = 'not a specific gene'
# 	discordances.loc[discordances['hugoSymbol'].isin(lung_genes), 'gvhd_filter'] = 'lung'
# 	discordances.loc[discordances['hugoSymbol'].isin(liver_genes), 'gvhd_filter'] = 'liver'
# 	discordances.loc[discordances['hugoSymbol'].isin(skin_genes), 'gvhd_filter'] = 'skin'
# 	discordances.loc[discordances['hugoSymbol'].isin(colon_genes), 'gvhd_filter'] = 'colon'
# 	discordances.loc[discordances['hugoSymbol'].isin(lacrimal_genes), 'gvhd_filter'] = 'lacrimal'
# 	discordances.loc[discordances['hugoSymbol'].isin(oral_genes), 'gvhd_filter'] = 'oral'
# 	discordances.loc[discordances['hugoSymbol'].isin(hemato_genes), 'gvhd_filter'] = 'hemato'

# 	lung_expr = lung.loc[lung['genes'].isin(lung_specific_genes)]
# 	liver_expr = liver.loc[liver["genes"].isin(liver_specific_genes)]
# 	skin_expr = skin.loc[skin["genes"].isin(skin_specific_genes)]
# 	colon_expr = colon.loc[colon["genes"].isin(colon_specific_genes)]
# 	lacrimal_expr = lacrimal.loc[lacrimal["genes"].isin(lacrimal_specific_genes)]
# 	oral_expr = oral.loc[oral["genes"].isin(oral_specific_genes)]
# 	hemato_expr = hemato.loc[hemato["genes"].isin(hemato_specific_genes)]

# 	with pd.ExcelWriter("GvHD_SpecificGenes_Expression.xlsx") as writer:
# 		lung_expr.to_excel(writer, sheet_name="lung", index=False)
# 		liver_expr.to_excel(writer, sheet_name="liver", index=False)
# 		skin_expr.to_excel(writer, sheet_name="skin", index=False)
# 		colon_expr.to_excel(writer, sheet_name="colon", index=False)
# 		lacrimal_expr.to_excel(writer, sheet_name="lacrimal", index=False)
# 		oral_expr.to_excel(writer, sheet_name="oral", index=False)
# 		hemato_expr.to_excel(writer, sheet_name="hemato", index=False)
		
# 	return(discordances)
