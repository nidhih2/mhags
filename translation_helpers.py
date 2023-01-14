import pandas as pd 
pd.options.mode.chained_assignment = None
import re

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import helpers as helpers

colsp = ['hugoSymbol','chromosome','start','end','variantClassification','variantType','refAllele','annotationTranscript',
		'cDnaChange','codonChange','proteinChange']

def parse_cdna_change(data, column):
	outputs = [['mutstart','mutend','mutref','mutalt']]
	for index, row in data.iterrows():
		if row.variantType in ['SNP', 'ONP','TNP','DNP']:
			spl = re.split(r'\.|>', row.loc[column])
			ref = re.split(r'\d+', spl[1])[1]
			start = re.split(r'[A-Z]+', spl[1])[0]
			if len(ref) > 1:
				end = start+len(ref)-1
			else: end = start
			vec = [start, end, ref, spl[2]]
			if row.variantType=='TNP':
				print('TNP')
				#print(row.cDnaChange)
				#print(vec)
			outputs.append(vec)
		elif (row.variantType=='INS') & ('_' in row.loc[column]):
			spl = re.split(r'\.|ins|_', row.loc[column])
			vec = [spl[1],spl[2],'',spl[3]]
			outputs.append(vec)
		elif row.variantType=='INS':
			spl = re.split(r'\.|ins', row.loc[column])
			vec = [spl[1],spl[1],'',spl[2]]
			outputs.append(vec)
		elif (row.variantType=='DEL') & ('_' in row.loc[column]):
			spl = re.split(r'\.|del|_', row.loc[column])
			vec = [spl[1],spl[2],spl[3],'']
			outputs.append(vec)
		elif row.variantType=='DEL':
			spl = re.split(r'\.|del', row.loc[column])
			vec = [spl[1],spl[1],spl[2],'']
			outputs.append(vec)
		else:
			print('Something bad happened')
			print(row)
			assert(False)

	parsed = pd.DataFrame(outputs[1:], columns=outputs[0])
	return(parsed)

def parse_protein_change(data):
	outputs = [['protstart','protend','protref','protalt','fs']]

	for index, row in data.iterrows():
		fs = 'no'
		if row.variantType in['SNP', 'ONP','DNP','TNP']: ## TODO: DOUBLE CHECK THE ONP
			spl = re.split(r'\.|[0-9]+', row.proteinChange)
			spl2 = re.split(r'[A-Z]|\*', row.proteinChange)
			outputs.append([int(spl2[1]),int(spl2[1]),spl[1],spl[2],fs])

		elif row.variantType=='INS':
			if ('fs' in row.proteinChange) | ('SVP' in row.proteinChange):
				spl = re.split(r'\.|[0-9]+', row.proteinChange)
				spl2 = re.split(r'[A-Z]+|fs|\*', row.proteinChange)
				if 'fs' in row.proteinChange: 
					fs = 'yes'
					outputs.append([int(spl2[1]),int(spl2[1]),spl[1], '', fs])
				else:
					outputs.append([int(spl2[1]),int(spl2[1]), '', spl[1], fs])

			else:
				if '_' in row.proteinChange:
					spl = re.split(r'\.|_|ins', row.proteinChange)
					outputs.append([int(spl[1]),int(spl[2]),'',spl[3],fs])
				else: ### TODD: I NEED TO DEBUG THIS TO MAKE SURE IT WORKS PROPERLY -- AND DOUBLE CHECK WHAT THE LACK OF UNDERSCORE IS SUPPOSED TO PRODUCE
					spl = re.split(r'\.|[0-9]+', row.proteinChange)
					spl2 = re.split(r'[A-Z]+', row.proteinChange)

					outputs.append([int(spl2[1]), int(spl2[1]), spl[1], spl[2],'no'])


		elif row.variantType=='DEL':
			if '>' in row.proteinChange:
				spl = re.split(r'[0-9]+|>', row.proteinChange)
				spl2 = re.split(r'\.|[A-Z]+|_', row.proteinChange)
				outputs.append([int(spl2[1]),int(spl2[1]),spl[2],spl[3], fs])

			else:
				if ('*' in row.proteinChange) & ('fs' in row.proteinChange): ## @TODO: NEED TO DEBUG THIS TO MAKE SURE IT DOES WHAT I THINK IT DOES
					# print('frameshift deletion with start')
					# print(row.proteinChange)
					spl = re.split(r'\*|fs', row.proteinChange)
					# print(spl)
					outputs.append([int(spl[1]), int(spl[1]), '*','*', 'yes'])


				elif '*' in row.proteinChange:
					print('deletion with star case')
					print(row.proteinChange)
					print(row.cDnaChange)
					print(row.variantType)
					print(row.variantClassification)
					spl = re.split(r'\*|\.|[0-9]+|del|fs', row.proteinChange)
					spl2 = re.split(r'[A-Z]+|del|fs', row.proteinChange)
					print(spl)
					print(spl2)
					### FOR WEIRD EDGE CASE LOOKING LIKE: p.SHGAR*RS328del (deletion spanning stop codon)
					if len(spl2) > 3: print('deletion with *, continuing'); continue
					print()
					outputs.append([int(spl2[2]), int(spl2[2]), spl[1], spl[2],'no'])
					# assert(False)

				## this is the case for normal deletions
				else:
					# print('normal deletion')
					spl = re.split(r'\.|[0-9]+|del|fs', row.proteinChange)
					spl2 = re.split(r'[A-Z]+|del|fs', row.proteinChange)
					# print(row.proteinChange)
					# print(spl)
					# print(spl2)
					if 'fs' in row.proteinChange: fs = 'yes'
					outputs.append([int(spl2[1]),int(spl2[1]),spl[1],'', fs])

		else:
			print('neither option')

	parsed = pd.DataFrame(outputs[1:], columns=outputs[0])

	return(parsed)

def determine_mutation_order(data, sample):
	"""
	need to split the data table into two, or find a way to encode mutations that are in the same place
	1: deal with frameshifts that truncate the sequence on one strand, when there are still mutations for the other strand
	2: two mutations happening in the same codon that cannot be on the same strand

	"""
	prev_protindex = -1
	prev_cdnaindex = -1
	mutation_index = -1
	prev_type = ''
	prev_class = ''
	mutation_indices = []
	which_df = []
	## this is incremented every time there is a new mutation
	## sometimes there are two diferent mutations in the same place, in those cases it does not get incremented
	df1 = [list(data.columns)+[sample+'_mut_index']]
	df2 = [list(data.columns)+[sample+'_mut_index']]
	
	for index, row in data.iterrows():

		# print(row.cDnaChange)
		# print(row.proteinChange)
		# print(row.mutstart)
		# print(prev_cdnaindex)
		
		if row[sample+'_loci']=='hom':
			## if it's a double mutation, and it's not the same base, but both are snps, we want to keep both as usual
			if (row.protstart == prev_protindex) & (row.mutstart != prev_cdnaindex) & (row.variantType=='SNP') & (prev_type=='SNP'):
				# print('hom, two snps, adjacent mutations')
				# print('index: %i' %index)
				mutation_index += 1
				df1.append(list(row)+[mutation_index])
				df2.append(list(row)+[mutation_index])
				which_df.append(3)

			## two SNPs in the exact same place. for this case remove one from the one data frame and keep from the other (have to check that the previous was also homozygous), actually, this would (should) only happen in that case
			elif (row.protstart==prev_protindex) & (row.mutstart == prev_cdnaindex) & (row.variantType=='SNP') & (prev_type=='SNP'):
				print('\n\nTWO HOM SNPS IN THE EXACT SAME PLACE\n')
				## putting the mutation in the second dataframe
				print(data[colsp+[sample+'_GT',sample]])
				print(row[colsp+[sample+'_GT',sample]])
				df2 = df2[:len(df2)-1]
				df2.append(list(row)+[mutation_index])
				which_df.append(2)

			elif (row.protstart==prev_protindex) & (row.mutstart == prev_cdnaindex):
				print('\n\nTWO HOM MUTATIONS IN THE EXACT SAME PLACE\n')
				## putting the mutation in the second dataframe
				print(data[colsp+[sample+'_GT',sample]])
				print(row[colsp+[sample+'_GT',sample]])
				print('\nprevs')
				print(prev_type)
				print(prev_class)
				print(prev_protindex)
				print(prev_cdnaindex)
				# assert(False)
				df2 = df2[:len(df2)-1]
				df2.append(list(row)+[mutation_index])
				which_df.append(2)
			## generic case for the start location being the same as the previous, this accounts for cases with framshifts, these sould still be added in both locations
			## just add to the second dataframe
			## CHANGED: add to both first and second
			elif (row.protstart == prev_protindex) & (row.mutstart > prev_cdnaindex):
				mutation_index += 1
				df2.append(list(row)+[mutation_index])
				df1.append(list(row)+[mutation_index])
				which_df.append(3)

			## two mutations in the exact same place... this is probably wrong. check the quality score on them
			## this happens with SNPs sometimes because it will just compare reads of the referencd with reads of the alternate. so it's actually a het site with two different alleles
			

			# just for frameshifts,currently just adding to the second table, but I guess I should add to the first as well?
			## CHANGED: add to both first and second tables
			elif 'fs' in row.proteinChange:
				# print('\nhom, just a frameshift, stil putting in both tables\n')
				mutation_index += 1

				df2.append(list(row)+[mutation_index])
				df1.append(list(row)+[mutation_index])
				which_df.append(3)

				tmpdf = pd.DataFrame(df2[1:], columns=df2[0])

			else:
				mutation_index += 1
				df1.append(list(row)+[mutation_index])
				df2.append(list(row)+[mutation_index])
				which_df.append(3)

		elif row[sample+'_loci']=='het':

			## the same codon is changed, but not the same base, and both are SNPs keep together
			if (row.protstart == prev_protindex) & (row.mutstart != prev_cdnaindex) & (row.variantType=='SNP') & (prev_type=='SNP'):
				# print('het, two snps, adjacent mutations')
				# print('index: %i' %index)
				mutation_index += 1
				df1.append(list(row) + [mutation_index])
				which_df.append(1)

			## both SNPS and exact same base is changed, put on opposite transcripts, don't increment the counter because same mutation location
			elif (row.protstart == prev_protindex) & (row.mutstart == prev_cdnaindex) & (row.variantType=='SNP') & (prev_type=='SNP'):
				# print('het, two snps, exact same location')
				# print('index: %i' %index)
				df2.append(list(row) + [mutation_index])
				which_df.append(2)

			## for framsehift that happens in the same codon as a previous mutation
			elif ('fs' in row.proteinChange) & (row.protstart == prev_protindex):
				if (prev_class in ['FRAME_SHIFT_INS','FRAME_SHIFT_DEL']): ## this is rare but it does seem to happen
					print('\nTWO FRAME SHIFTS IN THE SAME PLACE!')
					df1.append(list(row)+[mutation_index])
					which_df.append(1)
				else: 
					# print('het, frameshift in same place as previous mutation')
					# print('index: %i' %index)
					df2.append(list(row)+[mutation_index])
					which_df.append(2)

			## IF AN INDEL HAPPENS IMMEDIATELY AFTER A SNP, THEY WILL BE PUT ON DIFFERENT STRANDS


			## for a normal frameshift, just put it on the second so it doesn't interfere with the bulk of the mutations
			### (should I just not add any more mutations to the end of this one then? will that just cause problems?)
			### (answer: yes. everything after a frameshift will be seen as homozygous because the frameshift messes with it. How do I deal with this?)
			### also use this for the case of NONSENSE MUTATIONS
			### 	since these cause a premature stop, they should't be on the same strand as the rest of the mutations	
			elif ('fs' in row.proteinChange) | (row.variantClassification == 'NONSENSE'):
				# print('het, frameshift or a nonsense')
				# print('index: %i' %index)
				## for the cases where there is another mutation in the same location in the previous spot, don't increment the counter
				if row.protstart != prev_protindex:
					mutation_index += 1
				df2.append(list(row)+[mutation_index])
				which_df.append(2)

			## master condition for a mutation happening in the same place as the previous one, when both are not SNPs
			elif row.protstart == prev_protindex: 
				# print('het, any other double mutation')
				df2.append(list(row)+[mutation_index])
				which_df.append(2)
				# assert(False)
			else:							## for all other het mutations add to the first sequence
				mutation_index += 1
				df1.append(list(row)+[mutation_index])
				which_df.append(1)
		else:
			print('something went wrong')
			print(data[colsp+[sample+'_GT']])
			print(row[sample+'_loci'])
			assert(False)
		prev_protindex = row.protstart
		prev_cdnaindex = row.mutstart
		prev_type = row.variantType
		prev_class = row.variantClassification
		mutation_indices.append(mutation_index)

	df1 = pd.DataFrame(df1[1:], columns=df1[0])
	df2 = pd.DataFrame(df2[1:], columns=df2[0])
	data[sample+'_mut_index'] = mutation_indices
	data[sample+'_whichdf'] = which_df
	# print(df1)
	# print(df2)

	return data, df1, df2


def find_mutation_indices(table, sample, num):

	"""
	Function for determining the location of the mutation in the mutated protein.
	This is necessary because in the case of an indel, the location provided in the hgvs might be off

	@param table: the data
	@param sample: host or donor
	@param num: whether it's the first or second protein

	@returns: a table with the 'mutation base index' for each of the mutations (the number /shift in index to add to the hgvs position)

	"""

	mutation_index_base = 0
	base_indices = []
	# base_indices = [mutation_index_base]

	if table.empty:
		return(pd.DataFrame(columns=['annotationTranscript','cDnaChange','proteinChange',sample+'_base_index_'+num]))
	for index, row in table.iterrows():
		# base_indices.append(mutation_index_base)

		if row.variantType == 'INS':
			mutation_index_base += len(row.protalt)
		elif row.variantType == 'DEL':
			mutation_index_base -= len(row.protref)

		base_indices.append(mutation_index_base)
		
	table[sample+'_base_index_'+num] = base_indices
	table = table[['annotationTranscript','cDnaChange','proteinChange',sample+'_base_index_'+num]]
	return(table)

def chop_sequence(data):

	segments = []
	startindex = 0
	seqpad = []

	wt_seq = data.loc[0,'mrna']
	for index, row in data.iterrows():
		# print('row %i' %index)
		if row.variantType in ['SNP','DNP','TNP']: ### @TODO: CHECK THAT THIS IS CORRECT OR ONP
			segment = wt_seq[startindex:int(row.mutstart)-1]
			startindex = int(row.mutend)
			segments.append(segment)
			try:
				assert(wt_seq[int(row.mutstart)-1]==row.mutref)
			except:
				continue

		elif row.variantType == 'INS':
			segment = wt_seq[startindex:int(row.mutstart)]
			startindex = int(row.mutstart)
			segments.append(segment)
		elif row.variantType == 'DEL':
			segment = wt_seq[startindex:int(row.mutstart)-1]
			startindex = int(row.mutend)
			segments.append(segment)
			try:
				assert(wt_seq[int(row.mutstart)-1:int(row.mutend)]==row.mutref)
			except:
				continue
		else:
			print('isnt a SNP, INS, or DEL')
			print(row[colsp])
			print(row.variantType)
			assert(False)

	segments.append(wt_seq[startindex:])

	return(segments)

def fetch_mutant_nucs(data):

	mut1 = []
	refs = []
	for index, row in data.iterrows():
		mut1.append(row.mutalt)

		refs.append(row.mutref)

	return(mut1, refs)

def apply_mutations(segments, mutations):
	# print(segments)
	# print(mutations)

	mutseq = ''
	mutseq_print = ''
	try:
		assert(len(segments) == len(mutations)+1)
	except:
		print('the number of segments and mutations don\'t align')
		print(segments)
		print(mutations)
		assert(False)

	for i in range(len(mutations)):
		mutseq += segments[i]
		mutseq += mutations[i]

		mutseq_print += segments[i]
		mutseq_print += ' '
		mutseq_print += mutations[i]
		mutseq_print += ' '

	mutseq += segments[-1]
	mutseq_print += segments[-1]

	return(mutseq)


def extract_full_proteins(data, table1, table2, gene, sample):
	genes_version = {'annotationTranscript':[], 
					sample+'_cdnamut1':[], 
					sample+'_cdnamut2':[], 
					sample+'_protmut1':[],
					sample+'_protmut2':[]}
	
	## deal with table 1
	## if the table is empty, then there was just a frameshift (probably) and it was put in table 2
	if table1.empty:
		## the rna and protein should be the reference 
		## because there are no rows, can't find the reference sequence. will be added in later
		mutated_sequences1, mrna1 = data.loc[0,'mrna'], data.loc[0,'mrna']
		chopped1 = [data.loc[0,'mrna']]
		mutated_prot1 = data.loc[0,'refProtein']
		mutations1, ref_alleles1 = [], []
		table1 = pd.DataFrame(columns=['hugoSymbol','annotationTranscript',sample+'_GT','proteinChange',sample+'_mut_index',sample+'_peptide1',sample+'_ref_peptide1'])
	else:
		chopped1 = chop_sequence(table1)
		mutations1, ref_alleles1 = fetch_mutant_nucs(table1)
		mutated_sequences1 = apply_mutations(chopped1, mutations1)
		mutated_prot1, mrna1 = helpers.translate(mutated_sequences1)

	if table2.empty:
		mutated_sequences2, mrna2 = data.loc[0,'mrna'], data.loc[0, 'mrna']
		chopped2 = [data.loc[0, 'mrna']]
		mutated_prot2 = data.loc[0,'refProtein']
		mutations2, ref_alleles2 = [], []
		table2 = pd.DataFrame(columns=['hugoSymbol','annotationTranscript',sample+'_GT','proteinChange',sample+'_mut_index',sample+'_peptide2',sample+'_ref_peptide2'])			
	else:
		chopped2 = chop_sequence(table2)
		mutations2, ref_alleles2 = fetch_mutant_nucs(table2)
		mutated_sequences2 = apply_mutations(chopped2, mutations2)
		mutated_prot2, mrna2 = helpers.translate(mutated_sequences2)

	genes_version['annotationTranscript'].append(gene)
	genes_version[sample+'_cdnamut1'].append(mrna1)
	genes_version[sample+'_cdnamut2'].append(mrna2)
	genes_version[sample+'_protmut1'].append(mutated_prot1)
	genes_version[sample+'_protmut2'].append(mutated_prot2)

	genes_version = pd.DataFrame(genes_version)
	return genes_version

def fetch_seqpad(data, sample, seqpad):

	print('\nFetching mutant peptides:')
	print(data)

	mutpeps1 = []
	mutpeps2 = []
	refpeps = []

	muts = []
	transcript = []

	# prev_base1 = 0
	# prev_base2 = 0
	base1 = 0
	base2 = 0
	next_base1 = 0
	next_base2 = 0

	for index, row in data.iterrows():
		wt_seq = row.refProtein

		prot1 = row[sample+'_protmut1']
		prot2 = row[sample+'_protmut2']

		if not pd.isnull(row[sample+'_base_index_1']):
			next_base1 = int(row[sample+'_base_index_1'])
		if not pd.isnull(row[sample+'_base_index_2']):
			next_base2 = int(row[sample+'_base_index_2'])


		if row.protstart > len(wt_seq):

			print(' ROW.PROTSTART > LEN(PROT1)')
			mutpep1 = '*'
		elif row.protstart > len(wt_seq):
			mutpep2 = '*'

		elif row.variantType in ['SNP','DNP','TNP']:


			start = max(0, row.protstart-seqpad)
			end = row.protend + seqpad

			refpeptide = wt_seq[start:end]

			mutation_index1 = base1+row.protstart-1
			mutation_index2 = base2+row.protstart-1

			if mutation_index1 > len(prot1):
				mutpep1 = '*'
			else:
				start1 = max(0, mutation_index1-seqpad)
				end1 = base1+end
				mutpep1 = prot1[start1:end1]

			if mutation_index2 > len(prot2):
				mutpep2 = '*'
			else:
				start2 = max(0, mutation_index2-seqpad)
				end2 = base2+end
				mutpep2 = prot2[start2:end2]

		elif row.variantType=='INS':

			value_add = 0
			if row.mutstart % 3 != 0: 
				value_add = int(len(row.mutalt)/3)

			start = max(0, row.protstart-seqpad)
			end = row.protend+seqpad+value_add
			refpeptide = wt_seq[start:end]

			mutation_index1 = base1+row.protstart
			start1 = max(0, mutation_index1-seqpad)

			mutation_index2 = base2+row.protstart
			start2 = max(0, mutation_index2-seqpad)

			if row[sample+'_whichdf'] in [1,3]:
				if row.variantClassification=='FRAME_SHIFT_INS':
					end1 = len(prot1) 
					start1 -= 1
				else:
					# print(row[['hugoSymbol','variantClassification','cDnaChange','proteinChange','mutstart','mutend','mutref','mutalt','protstart','protend','protref','protalt']])
					end1 = mutation_index1 + len(row.protalt) + seqpad + value_add
			else:
				# print('\nHERE')
				if row.variantClassification=='FRAME_SHIFT_INS':
					end1 = len(prot1) 
					start1 -= 1
				else:
					end1 = base1+end-value_add-1

			if row[sample+'_whichdf'] in [2,3]:
				if row.variantClassification=='FRAME_SHIFT_INS':
					end2 = len(prot2)
					start2 -= 1
				else:
					end2 = mutation_index2 + len(row.protalt) + seqpad + value_add
			else:
				if row.variantClassification=='FRAME_SHIFT_INS':
					end2 = len(prot2) 
					start2 -= 1
				else:
					end2 = base2+end-value_add-1

			mutpep1 = prot1[start1:end1]
			mutpep2 = prot2[start2:end2]


		elif row.variantType=='DEL':

			value_add = 0
			if (row.mutstart-1) % 3 != 0: 
				value_add = int(len(row.mutref)/3)-1

			start = max(0, row.protstart-seqpad)
			end = row.protend+seqpad+value_add
			refpeptide = wt_seq[start:end]

			if row[sample+'_whichdf'] in [1,3]:
				mutation_index1 = base1+row.protstart-1
				start1 = max(0, mutation_index1-seqpad)
				end1 = mutation_index1 + seqpad + value_add
			else:
				mutation_index1 = base1+row.protstart-1
				start1 = max(0, mutation_index1-seqpad)
				end1 = mutation_index1 + seqpad + value_add + len(row.protref)

			if row[sample+'_whichdf'] in [2,3]:
				mutation_index2 = base2+row.protstart-1
				start2 = max(0, mutation_index2-seqpad)
				end2 = mutation_index2 + seqpad + value_add
			else:
				mutation_index2 = base2+row.protstart-1
				start2 = max(0, mutation_index2-seqpad)
				end2 = mutation_index2 + seqpad + value_add + len(row.protref)

			if row.variantClassification == 'FRAME_SHIFT_DEL':
				end1 = len(prot1)
				end2 = len(prot2)

			mutpep1 = prot1[start1:end1]
			mutpep2 = prot2[start2:end2]


		else:
			print('mutation neither INS, DEL, nor SNP')
			print(row[['hugoSymbol','annotationTranscript','variantType','variantClassification','proteinChange','cDnaChange']])

			assert(False)

		if (row.host_hetsite == 'yes') & (row.host_whichdf in [1,3]):
			mutpep2 = '*'

		if (row.host_hetsite == 'yes') & (row.host_whichdf in [2,3]):
			mutpep1 = '*'

		base1 = next_base1
		base2 = next_base2
				
		muts.append(row.proteinChange)
		transcript.append(row.annotationTranscript)
		mutpeps1.append(mutpep1)
		mutpeps2.append(mutpep2)
		refpeps.append(refpeptide)

	peps = pd.DataFrame.from_dict({'annotationTranscript':transcript,
								'proteinChange':muts,
								sample+'_seqpad1':mutpeps1,
								sample+'_seqpad2':mutpeps2,
								sample+'_ref_peptide':refpeps})

	return(peps)