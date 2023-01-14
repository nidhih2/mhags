import pandas as pd 
import argparse

# Donor = the specific_discordances = tissue_specific_discordances.txt -> result = donor_discordances_after_blast.txt
# Host = the specific_discordances = donor_discordances_after_blast.txt -> result = host_discordances_after_blast.txt

# For donor -> no genes to remove (all are fair game)
# For host -> remove donor genes 
def filter_blast_search(specific_discordances, blast, sample_type, sample, condition, genes_to_remove=[]):
    specific_discordances = pd.read_csv(specific_discordances, sep="\t")
    blast = blast[(blast.peptide_len==blast.alignment_length) & (blast.p_ident==100)]
	## or the host situation, remove the exact matches in genes we care about. Just look at exact matches from other genes 
	## (these are genes likey expressed in other tissues, so we want to remove those, since they might be involved in GvHD)
    blast = blast[~blast.protein_gene.isin(genes_to_remove)]
    discordances = specific_discordances[~specific_discordances.pep.isin(blast.peptide_seq)]
    discordances.to_csv(sample + "_discordances_after_blast.txt", sep="\t")
    return(discordances)

def main():
    # Adding column name to the blast_results (sample_blastp_peptides_out.csv) 
    blast_results = pd.read_csv(blastp_result, sep=",", header=None) # result of blastp 
    blast_results.columns=['protein_id','peptide_id','protein_len','peptide_len','protein_start','protein_end','peptide_start','peptide_end','protein_seq','peptide_seq','evalue','score','alignment_length','p_ident','n_ident']
    blast_results = blast_results.assign(protein_gene=blast_results.protein_id.str.split('|', expand=True)[0])
    blast_results = blast_results.assign(protein_trans=blast_results.protein_id.str.split('|', expand=True)[2]) # protein transcript
    blast_results = blast_results.assign(peptide_gene=blast_results.peptide_id.str.split('|', expand=True)[0])
    blast_results = blast_results.assign(peptide_trans=blast_results.peptide_id.str.split('|', expand=True)[2])

    # specific_discordances = tissue_specific_discordances.txt
    if sample_type == "donor":
        donor_discordances = filter_blast_search(specific_discordances=specific_discordances, blast=blast_results, genes_to_remove=[], sample_type = "donor", sample = sample, condition=condition)
   
    # specific_discordances = donor_discordances_after_blast.txt with genes to remove 
    if sample_type == "host":
        donor_discordances = filter_blast_search(specific_discordances=specific_discordances, blast=blast_results, genes_to_remove=[], sample_type = "donor", sample = sample, condition=condition)        
        exclude_genes = set(donor_discordances.hugoSymbol)
        host_discordances = filter_blast_search(specific_discordances=specific_discordances, blast=blast_results, genes_to_remove=exclude_genes, sample_type = "host", sample = sample, condition=condition)
        host_discordances.to_csv('allgenes_pumas.txt', sep='\t', index=False)

        if condition == "GvL":
            final_dict = {"post_host_BLAST_GvL_peptides":len(set(host_discordances.pep)), "post_donor_BLAST_GvL_peptides":len(set(donor_discordances.pep)),
            "post_BLAST_GvL_genes":len(set(host_discordances["hugoSymbol"]))}
            final_dict_df = pd.DataFrame(final_dict, index=[0])
            final_dict_df.to_csv(condition+"_postBlast_dict_nums.csv", sep=",", index=False) # GvL_dict_nums.csv, GvHD_dict_nums.csv
        else:
            final_dict = {"post_host_BLAST_GvHD_peptides":len(set(host_discordances.pep)), "post_donor_BLAST_GvHD_genes":len(set(donor_discordances.pep))}
            final_dict_df = pd.DataFrame(final_dict, index=[0])
            final_dict_df.to_csv(condition+"_postBlast_dict_nums.csv", sep=",", index=False) # GvL_dict_nums.csv, GvHD_dict_nums.csv


if __name__ == "__main__":
    mutcount_cols = ['hugoSymbol','CHROM','start','end','refAllele','tumorSeqAllele1','tumorSeqAllele2']
    numbers_dict = {}

    parser = argparse.ArgumentParser()

    parser.add_argument("-specific_discordances", help='Tissue specific discordances to create peptides files/chr for blast search',default=None)
    parser.add_argument("-blastp_result", help="Peptides that are matched against the sample database", default=None)
    parser.add_argument("-sample_type", help = "Differentiates between donor and host")
    parser.add_argument("-sample", help = "The case name of donor and host for the resulting file")
    parser.add_argument("-condition", help="GvL or GvHD")
    args = parser.parse_args()
    
    specific_discordances = args.specific_discordances
    blastp_result = args.blastp_result
    sample_type = args.sample_type
    sample = args.sample
    condition = args.condition

    main()