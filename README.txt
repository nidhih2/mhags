Folder: mhags-pipeline 

####### Content:
1. 8 python scripts + 1 test.py 
2. 1 DockerFile
3. 1 req.txt 

####### Description: 
The pipeline takes as input Donor and Patient WES samples, Sex, Patient HLA typing (tested on HLA-Matched Donor and Patients) - performs Variant Calling using DeepVariant, Variant Annotation using Funcotator, Merges Donor and Patient annotated data simultaneous *construction of Donor and Patient Custom Proteome, Gene/Variant Filtration steps: 1) Retains Variants specific to Patients, 2) Looks for Non-Synonymous Variants (Protein AA seq altering variants - Missense, Non-stop, Nonsense, Frame Shifts), 3) Looks for Variants present in "Tissues of Interest" (Lung, Liver, Skin, Colon, Lacrimal and Oral Mucosa - **scRNA Seq exp). ***Peptides - checks if remnant peptides are present in Donor or Patient Custom Proteomes, if not present then predicts binding stability using HLAthena with the Patient HLAs (Peptide-HLA binding prediction), HLAthena post processing: picks out peptides that have a binding stability of <0.5

#######
*construction of Donor and Patient Custom Proteome = Explained in presentation

**scRNA Seq exp = For every tissue (organ) = 
1. Removed immune cells 
2. Reclustered and annotated tissue specific cell types 
3. Strategically eliminated drop-out effect 
(Path to scRNA Exp: mhags-pipeline-work/sc_analysis) 

***Peptides = Peptides are taken from Funcotator Transcript subfolder (after Variant annotation)

####### Storage details and other info:
1. Docker name: "nidhihookeri/minors-pipeline-2" (docker.com)
2. TERRA Workspace name: broad-fireclous-wuclonal/mHAgs_pipeline_nidhi 
3. TERRA Workflow name/Firecloud: mhags_pipeline 
4. GCP Bucket info: wu-lab-archives/wld-instance-2/mhags-data (gcloud compute ssh wld-instance-2 --project wu-lab-archives --zone us-central1-a) 
5. Presentation: mhags-pipeline-work/mhags_pipeline.pptx
6. Final WDL: mhags-pipeline-work/WDL/mhags_pipeline.15.wdl
7. Associated inputs JSON:  mhags-pipeline-work/WDL/mhags_pipeline.15.json

####### ORDER OF RUN (IF LOCALLY): 
1. bmt-simulation.py (includes as child files: functions.py, helpers.py, translation_helpers.py)
2. blast_preprocessing.py
3. blast_postprocessing.py
4. HLAthena_preprocessing.py
5. HLAthena_postprocessing.py

####### USAGE: 

1. 
python3 mhags-pipeline/bmt-simulation.py -donor_vcf /path/to/donorVCF/from/Funcotator -host_vcf /path/to/hostVCF/from/Funcotator -host_sex M/F -donor_sex M/F -tissue AML/CML/CLL -gvhd yes/no 

OUTPUT  CASES TAKEN INTO CONSIDERATION = (M->M, M->F, F->F): 
A. AML_specific_discordances.txt
B. GvHD_specific_discordances.txt

2. 
2A. python3 mhags-pipeline/blast_preprocessing.py -specific_discordances /for/GvL/AML_specific_discordances.txt 
OUTPUT: 
A. peptides_for_blast.fa (FASTA file for GvL output)

2B. python3 mhags-pipeline/blast_preprocessing.py -specific_discordances /for/GvL/GvHD_specific_discordances.txt 
OUTPUT: 
A. peptides_for_blast.fa (FASTA file for GvHD output)

3. (Maintain sequence of Donor to Host - important)
# DONOR GvL
3Aa.  python3 mhags-pipeline/blast_postprocessing.py -specific_discordances /for/GvL/AML_specific_discordances.txt -blastp_result /after/makeblastdb/and/blastp/blastp_result.csv -sample donor_basename -sample_type "donor" -condition GvL

3Ab MAKEBLASTDB AND BLASTP

# HOST GvL
3Ac.  python3 mhags-pipeline/blast_postprocessing.py -specific_discordances /discordances/from/3Aa -blastp_result /after/makeblastdb/and/blastp/after/donor/blastp_result.csv -sample host_basename -sample_type "host" -condition GvL

# DONOR GvHD
3Ba.  python3 mhags-pipeline/blast_postprocessing.py -specific_discordances /for/GvHD/GvHD_specific_discordances.txt -blastp_result /after/makeblastdb/and/blastp/blastp_result.csv -sample donor_basename -sample_type "donor" -condition GvL

3Bb MAKEBLASTDB AND BLASTP

# HOST GvHD
3Bc  python3 mhags-pipeline/blast_postprocessing.py -specific_discordances /discordances/from/3Ba -blastp_result /after/makeblastdb/and/blastp/after/donor/blastp_result.csv -sample host_basename -sample_type "host" -condition GvL

4. 
4A  python3 mhags-pipeline/HLAthena_preprocessing.py -discordances_after_blast /discordances/from/3Ac -autosomal "True"
4B  python3 mhags-pipeline/HLAthena_preprocessing.py -discordances_after_blast /discordances/from/3Bc -autosomal "True"

5. HLAthena 

6. 
6A  python3 mhags-pipeline/HLAthena_postprocessing.py -sample_predictions /HLAthena/output/sample_predictions -discordances_after_blast /discordances/output/from/3Ac/for/GvL

6B  python3 mhags-pipeline/HLAthena_postprocessing.py -sample_predictions /HLAthena/output/sample_predictions -discordances_after_blast /discordances/output/from/3Bc/for/GvHD

FINAL OUTPUTS: 
A. GvL HLAthena postprocessing binding putative minors antigens (text file)
B. GvHD HLAthena postprocessing binding putative minors antigens (text file)



