import pandas as pd
import argparse
# run 1 - specific_discordances as input 
# Takes the pep and the peptide_fasta_id from the specific_discordances and writes it into a fasta format 
#eg line 1: pep_fasta_id ; line 2 - pep (foreign peptide)
def blast_files(specific_discordances):
    discordances = pd.read_csv(specific_discordances, sep="\t")
    for chrom in pd.unique(discordances.CHROM):
        discs = discordances[discordances.CHROM==chrom]
        print(chrom)
        print('Number of peptides in %s dataset: %i' %(chrom, discs.shape[0]))
        discs = discs[['peptide_fasta_id','pep']].dropna().drop_duplicates('pep')
        print('Number of peptides after dropping duplicates and null values %s: %i' %(chrom, discs.shape[0]))

        print('Writing fa file of %s peptides for blast search...' %chrom, "\n")

        name = chrom
        with open('peptides_for_blast.fa','a') as f:
            for index, row in discs.iterrows():
                if len(row.pep) >= 8 & len(row.pep) <= 11:
                    f.write('>'+row.peptide_fasta_id)
                    f.write("\n")
                    f.write(row.pep)
                    f.write("\n")
        f.close()
    

def main():
    blast_files(specific_discordances = specific_discordances)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-specific_discordances", help='Tissue specific discrodances to create peptise files/chr for blast search',default=None)

    args = parser.parse_args()
    specific_discordances = args.specific_discordances
    main()



