from __future__ import annotations
from audioop import add
import dis
import pandas as pd
import argparse


def process_hlathena_predictions(hlathena_predictions, discordances, id_col, sample_name):

    hlathena_predictions = pd.read_csv(hlathena_predictions, sep="\t")
    discordances = pd.read_csv(discordances, sep="\t")

    hlathena_predictions[id_col] = hlathena_predictions[id_col].str.split(';')
    hlathena_predictions = hlathena_predictions.explode(id_col)

    ## first find reference predictions:
    mut_predictions = hlathena_predictions[hlathena_predictions.source=='host']
    ref_predictions = hlathena_predictions[hlathena_predictions.source=='donor']

    # taking only the important columns into consideration
    hlathena_mut_want = mut_predictions[['pep',id_col, 'len','assign.MSi_ranks', 'assign.MSi_allele']]
    ## take just the peptide and prediction columns from the reference predictions, to add in to the data object
    hlathena_ref_want = ref_predictions[['pep',id_col, 'best.MSi']]
    hlathena_ref_want.columns = ['ref_pep', id_col, 'ref_pep_best.MSi']

    # removing the unknown and any predictions have a threshold greater than 0.5 (strong binders)
    hlathena_mut_want = hlathena_mut_want[hlathena_mut_want.loc[:,'assign.MSi_allele']!='unknown']
    hlathena_mut_want = hlathena_mut_want[hlathena_mut_want.loc[:,'assign.MSi_ranks'] < 0.5]

    hlathena_weak_alldata = hlathena_ref_want.merge(hlathena_mut_want, on=[id_col], how='inner')
    discordances = discordances.merge(hlathena_weak_alldata, how='inner')	

    ## MERGE WITH THE APPROPRIATE REFERENCE PEPTIDE PREDICTION
    discordances = discordances.merge(hlathena_ref_want, how='left')

    return discordances

# Redundant function to retain the tissue of interest that is used in the workflow: used to extract the "tissue_type:AML/ALL/BLOOD etc" column from the discordances data
def tissue_of_interest(tissue_type):
    return tissue_type


def main():
    hla_binding = process_hlathena_predictions(hlathena_predictions=hlathena_predictions, discordances=discordances, id_col="transid", sample_name=sample_name)
    type_tissue = tissue_of_interest(tissue_type=tissue_type)

    cols_keep = ['hugoSymbol','CHROM','cDnaChange','variantType', 'proteinChange','pep_length', 'pep', 'ref_pep',
                 'assign.MSi_ranks', 'ref_pep_best.MSi', 'assign.MSi_allele','annotationTranscript']

    additional_cols = [tissue_type, "gvhd_classfication", "BLOOD"] # add patient_exp if needed in the future (comes along with GvL?)
    for col in additional_cols:
        if col in hla_binding.columns:
            cols_keep.append(col)

    hla_binding = hla_binding[cols_keep]

    # there are isoforms of the a feew genes, all the annotation Transcripts are strung into one line
    annotationTrans = pd.DataFrame(hla_binding.groupby(["hugoSymbol"]).agg(annotationTranscript=('annotationTranscript', ",".join)))
    hla_binding = hla_binding.drop_duplicates(subset="hugoSymbol")
    hla_binding.drop(labels="annotationTranscript", axis=1, inplace=True)
    merged = hla_binding.merge(annotationTrans, on="hugoSymbol", how="inner")

    # sorting the values wrt peptide length and the MSi_allele
    merged = merged.sort_values(by=["pep_length", "assign.MSi_allele"])

        # GvL log file entry 
    if tissue_type in hla_binding.columns:
        tissue_discs = hla_binding[hla_binding[tissue_type]=='yes'] # function tissue_type returns tissue used in the workflow as str 
        blood_discs = hla_binding[hla_binding.BLOOD=='yes']

        final_dict_gvl = {"GvL_predicted_binders_tot":len(set(merged["pep"])), 
        "GvL_predicted_binders_AML_filter":len(set(tissue_discs.pep)),
        "GvL_predicted_binders_heme_filter":len(set(blood_discs.pep))}
        final_dict_gvl_df = pd.DataFrame(final_dict_gvl, index=[0])
        final_dict_gvl_df.to_csv("GvL_postHLAthena_dict_nums.csv", sep=",", index=False)

        # GvHD log file entry
    if "gvhd_classfication" in hla_binding.columns:
        gvhd_discs_gene = hla_binding[["gvhd_classfication", "hugoSymbol"]].drop_duplicates()
        gvhd_discs = gvhd_discs_gene["gvhd_classfication"].value_counts().to_dict()

        gvhd_discs_gene["gvhd_classfication_list"] = gvhd_discs_gene["gvhd_classfication"].str.split(",")
        gvhd_tissue_bins = gvhd_discs_gene["gvhd_classfication_list"].str.len().value_counts().to_dict()

        final_dict_gvhd = {"GvHD_predicted_binders_tot":len(set(merged["pep"])),
        "GvHD_predicted_binders_pan_expr":gvhd_discs["pan expressed gene"], "GvHD_predicted_binders_non_heme_expr":gvhd_discs["lung,liver,skin,colon,lacrimal,oral"],
        "GvHD_predicted_binders_lung":gvhd_discs["lung"], "GvHD_predicted_binders_liver":gvhd_discs["liver"], "GvHD_predicted_binders_skin":gvhd_discs["skin"],
        "GvHD_predicted_binders_GI":gvhd_discs["colon"], "GvHD_predicted_binders_lacrimal":gvhd_discs["lacrimal"],"GvHD_predicted_binders_oral":gvhd_discs["oral"],
        "GvHD_predicted_binders_heme":gvhd_discs["hematopoietic"],	"GvHD_predicted_binders_2_tissues":gvhd_tissue_bins[2],	"GvHD_predicted_binders_3_tissues":gvhd_tissue_bins[3],
        "GvHD_predicted_binders_4_tissues":gvhd_tissue_bins[4],"GvHD_predicted_binders_5_tissues":gvhd_tissue_bins[5]}
        final_dict_gvhd_df = pd.DataFrame(final_dict_gvhd, index=[0])
        final_dict_gvhd_df.to_csv("GvHD_postHLAthena_dict_nums.txt", index=False)

    # Re ordering and renaming the files in the final output
    gvl_cols = ['gene', 'CHROM', 'DNA-Change', 'Variant-Type', 'Protein-Change', 'Peptide-Length', 'Variant-Peptide', 'Reference-Peptide',  
                         'Variant-Peptide-MSi-HLAthena-Rank', 'Reference-Peptide-Best-MSi', 'Variant-Peptide-MSi-HLAthena-Allele', 'Transcript-Annotation']
    gvhd_cols = ['gene', 'CHROM', 'DNA-Change', 'Variant-Type', 'Protein-Change', 'Peptide-Length', 'Variant-Peptide', 'Reference-Peptide',  
                         'Variant-Peptide-MSi-HLAthena-Rank', 'Reference-Peptide-Best-MSi', 'Variant-Peptide-MSi-HLAthena-Allele','GvHD-Tissues','Transcript-Annotation']

    # renaming the columns in the final files for GvL and GvHD
    if tissue_type in hla_binding.columns:
        cols_keep_gvl = ['hugoSymbol','CHROM','cDnaChange','variantType', 'proteinChange','pep_length', 'pep', 'ref_pep',  'assign.MSi_ranks', 'ref_pep_best.MSi', 'assign.MSi_allele', 'annotationTranscript']
        merged = merged[cols_keep_gvl]
        merged.columns = gvl_cols
    else: 
        cols_keep_gvhd = ['hugoSymbol','CHROM','cDnaChange','variantType', 'proteinChange','pep_length','pep', 'ref_pep',  'assign.MSi_ranks', 'ref_pep_best.MSi', 'assign.MSi_allele', 'gvhd_classfication', 'annotationTranscript']
        merged = merged[cols_keep_gvhd]
        merged.columns = gvhd_cols

    # final output file name to be processed and presented in TERRA
    file_name = sample_name + "_binding_putativeMinorAntigens.txt"
    merged.to_csv(file_name, sep='\t', index=False)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-hlathena_predictions', help = "sample_predictions from HLAthena output separate for GvL and GvHD")
    parser.add_argument("-discordances", help = "allgenes_puma.txt")
    parser.add_argument("-tissue_type", help = "Redundant way of calling the tissue of interest for each case") # USED FOR THE REDUNDANT FUNCTION: tissue_type
    parser.add_argument("-sample_name")

    args = parser.parse_args()
    hlathena_predictions = args.hlathena_predictions
    discordances = args.discordances
    tissue_type = args.tissue_type # USED FOR THE REDUNDANT FUNCTION: tissue_type
    sample_name = args.sample_name

    main()