import "https://api.firecloud.org/ga4gh/v1/tools/HLAthena_v1:HLAthena_v1_external/versions/8/plain-WDL/descriptor" as HLAthena

workflow realign {
  
  String sample_id
  String donor_output_basename
  String host_output_basename

  File? donor_input_cram
  File? host_input_cram

  File? donor_input_bam
  File? donor_input_bai
  File? host_input_bam
  File? host_input_bai

  # hg38 ref files
  File original_ref_fasta 
  File original_ref_fasta_index

  # hg19 ref files 
  File ref_fasta
  File ref_fasta_index

  # BWA files 
  File ref_dict
  File ref_bwt
  File ref_amb
  File ref_ann
  File ref_pac
  File ref_sa

  # Single Nucleotide Polymorphism Database - hg19
  File dbSNP_vcf
  File dbSNP_vcf_index

  # Deletion/insertion polymorphism VCF - hg19 
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices

  # Deepvariant input 
  File? capture_interval_list
  Int runtime_disk = 1
  String runtime_docker

  # Functotator input
  String? data_sources_tar_gz

  # bmt input 
  String tissue
  String donor_sex
  String host_sex
  String gvhd_filter

  # HLAthena input
  File hlatypes_file
  Boolean aggregate_pep
  String peptide_col_name
  Boolean exists_expr
  String expr_col_name
  Array[String] lens
  Boolean exists_ctex

  # Disks and multipliers
  Int? additional_disk
  Int preemptible_tries
  Int agg_preemptible_tries

  Float cram_disk_multiplier = 8
  Float bwa_disk_multiplier = 2.5
  Int compression_level = 2
  Float md_disk_multiplier = 2.25
  Float sort_sam_disk_multiplier = 10

  Float original_ref_size = size(original_ref_fasta, "GB") + size(original_ref_fasta_index, "GB")
  Int disk_pad = select_first([additional_disk, 20]) 

  # Output basenames for different tasks 
  #String? donor_output_basename = basename(donor_input_cram, ".cram")
  #String? host_output_basename = basename(host_input_cram, ".cram")

  String? donor_recalibrated_bam_basename = donor_output_basename + ".aligned.duplicates_marked.recalibrated"
  String? host_recalibrated_bam_basename = host_output_basename + ".aligned.duplicates_marked.recalibrated"


  # BWA tool string 
  # -p: prefix for the database (same as the file name); -v: verbose level of output; -t: threads; -k: min seed length 
  String bwa_commandline="bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"

  # File sizes 
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Float bwa_ref_size = ref_size + size(ref_amb, "GB") + size(ref_ann, "GB") + size(ref_bwt, "GB") + size(ref_pac, "GB") + size(ref_sa, "GB")
  Float dbsnp_size = size(dbSNP_vcf, "GB") + size(dbSNP_vcf_index, "GB")

  # CRAM to BAM conversion using samtools; slow due to intermediate file output of SAM file 

  if(defined(host_input_cram) &&  defined(donor_input_cram) && !defined(host_input_bam) && !defined(donor_input_bam)) {

    call CramToBam as donor_CramToBam {
      input:
        cram_file = donor_input_cram,
        output_basename = donor_output_basename,
        ref_fasta = original_ref_fasta,
        ref_fasta_index = original_ref_fasta_index,
        disk_size = ceil(size(donor_input_cram, "GB") * 6 + original_ref_size) + disk_pad
    }

    call CramToBam as host_CramToBam {
      input:
        cram_file = host_input_cram,
        output_basename = host_output_basename,
        ref_fasta = original_ref_fasta,
        ref_fasta_index = original_ref_fasta_index,
        disk_size = ceil(size(host_input_cram, "GB") * 6 + original_ref_size) + disk_pad
    }
  }

  File? donor_bam_data = select_first([donor_CramToBam.output_bam, donor_input_bam])
  File? donor_bai_data = select_first([donor_CramToBam.output_bam_index, donor_input_bai])
  File? host_bam_data = select_first([host_CramToBam.output_bam, host_input_bam])
  File? host_bai_data = select_first([host_CramToBam.output_bam_index, host_input_bai])

  # GATK tool; BAM to uBAM
  Float donor_input_size = size(donor_bam_data, "GB")
  call RevertSam as donor_RevertSam {
    input:
      input_bam = donor_bam_data,
      disk_size = ceil(donor_input_size * cram_disk_multiplier) + disk_pad,
      output_basename = donor_output_basename
  }

  Float host_input_size = size(host_bam_data, "GB")
  call RevertSam as host_RevertSam {
    input:
      input_bam = host_bam_data,
      disk_size = ceil(host_input_size * cram_disk_multiplier) + disk_pad,
      output_basename = host_output_basename
  }

  # GATK tool; sorts the uBAM wrt query name (important to group by queryname for MergeBamAlignment task - recommended by GATK)

  String donor_sortSam_output_basename = basename(donor_RevertSam.unmapped_bam, ".bam")
  Float donor_revertSam_bam_size = size(donor_RevertSam.unmapped_bam, "GB")
  call SortSam as donor_SortSam {
    input:
      input_bam = donor_RevertSam.unmapped_bam,
      sorted_bam_name = donor_sortSam_output_basename + ".unmapped.bam",
      disk_size = ceil(donor_revertSam_bam_size * (10 * 1.75)) + disk_pad
  }

  String host_sortSam_output_basename = basename(host_RevertSam.unmapped_bam, ".bam")
  Float host_revertSam_bam_size = size(host_RevertSam.unmapped_bam, "GB")
  call SortSam as host_SortSam {
    input:
      input_bam = host_RevertSam.unmapped_bam,
      sorted_bam_name = host_sortSam_output_basename + ".unmapped.bam",
      disk_size = ceil(host_revertSam_bam_size * (10 * 1.75)) + disk_pad
  }

  # 1. Sam to fastq (GATK) conversion, 2. BWA (third-party tool) of sorted uBAM (aligning to hg38), 
  # 3. MergeBamAlignment (GATK) - uBAM data (meta-data: sample alias, library, barcodes, etc) + BWA:BAM alignment data -> mapped into a SAM file -> BAM (queryname sorted)

  Float donor_unmapped_bam_size = size(donor_SortSam.sorted_bam, "GB")
  call SamToFastqAndBwaMemAndMba as donor_SamToFastqAndBwaMemAndMba {
    input:
      input_bam = donor_SortSam.sorted_bam,
      bwa_commandline = bwa_commandline,
      output_bam_basename = donor_output_basename + ".aligned.unsorted",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_bwt = ref_bwt,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_pac = ref_pac,
      ref_sa = ref_sa,
      bwa_version = "0.7.15-r1140",
      disk_size = ceil(donor_unmapped_bam_size + bwa_ref_size + (bwa_disk_multiplier * donor_unmapped_bam_size)) + disk_pad,
      compression_level = compression_level,
      preemptible_tries = preemptible_tries
  }

  Float host_unmapped_bam_size = size(host_SortSam.sorted_bam, "GB")
  call SamToFastqAndBwaMemAndMba as host_SamToFastqAndBwaMemAndMba {
    input:
      input_bam = host_SortSam.sorted_bam,
      bwa_commandline = bwa_commandline,
      output_bam_basename = host_output_basename + ".aligned.unsorted",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_bwt = ref_bwt,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_pac = ref_pac,
      ref_sa = ref_sa,
      bwa_version = "0.7.15-r1140",
      disk_size = ceil(host_unmapped_bam_size + bwa_ref_size + (bwa_disk_multiplier * host_unmapped_bam_size)) + disk_pad,
      compression_level = compression_level,
      preemptible_tries = preemptible_tries
  }

  # Marks duplicate reads (defined as originating from a single fragment of DNA)
  call MarkDuplicates as donor_MarkDuplicates {
    input:
      input_bams = donor_SamToFastqAndBwaMemAndMba.output_bam,
      output_bam_basename = donor_output_basename + ".aligned.unsorted.duplicates_marked",
      metrics_filename = donor_output_basename + ".duplicate_metrics",
      disk_size = ceil(md_disk_multiplier * 40) + disk_pad,
      compression_level = compression_level,
      preemptible_tries = agg_preemptible_tries
  }

  call MarkDuplicates as host_MarkDuplicates {
    input:
      input_bams = host_SamToFastqAndBwaMemAndMba.output_bam,
      output_bam_basename = host_output_basename + ".aligned.unsorted.duplicates_marked",
      metrics_filename = host_output_basename + ".duplicate_metrics",
      disk_size = ceil(md_disk_multiplier * 40) + disk_pad,
      compression_level = compression_level,
      preemptible_tries = agg_preemptible_tries
  }

  # Sort aggregated+deduped BAM file wrt co-ordinates and calculates NM, 
  # Sorting the MarkDuplicated.BAM file wrt co-ordinates 

  Float donor_agg_bam_size = size(donor_MarkDuplicates.output_bam, "GB")
  call SortSampleBam as donor_SortSampleBam {
    input:
      input_bam = donor_MarkDuplicates.output_bam,
      output_bam_basename = donor_output_basename + ".aligned.duplicate_marked.sorted",
      disk_size = ceil(sort_sam_disk_multiplier * donor_agg_bam_size) + disk_pad,
      compression_level = compression_level,
      preemptible_tries = agg_preemptible_tries
  }

  Float host_agg_bam_size = size(host_MarkDuplicates.output_bam, "GB")
  call SortSampleBam as host_SortSampleBam {
    input:
      input_bam = host_MarkDuplicates.output_bam,
      output_bam_basename = host_output_basename + ".aligned.duplicate_marked.sorted",
      disk_size = ceil(sort_sam_disk_multiplier * host_agg_bam_size) + disk_pad,
      compression_level = compression_level,
      preemptible_tries = agg_preemptible_tries
  }

  # Create genomic intervals 
  # created intervals wrt sequence length
  call CreateSequenceGroupingTSV {
    input:
      ref_dict = ref_dict,
      preemptible_tries = preemptible_tries
  }

  # Perform BaseQuality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call BaseRecalibrator as donor_BaseRecalibrator {
      input:
        input_bam = donor_SortSampleBam.output_bam,
        input_bam_index = donor_SortSampleBam.output_bam_index,
        recalibration_report_filename = donor_output_basename + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        disk_size = ceil(donor_agg_bam_size + ref_size + dbsnp_size) + disk_pad,
        preemptible_tries = agg_preemptible_tries
    }
  }

  # Merge the recalibration reports resulting from by-interval recalibration
  # The reports are always the same size
  call GatherBqsrReports as donor_GatherBqsrReports {
    input:
      input_bqsr_reports = donor_BaseRecalibrator.recalibration_report,
      output_report_filename = donor_output_basename + ".recal_data.csv",
      disk_size = disk_pad,
      preemptible_tries = preemptible_tries
  }

  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {
    # Apply the recalibration model by interval
    call ApplyBQSR as donor_ApplyBQSR {
      input:
        input_bam = donor_SortSampleBam.output_bam,
        input_bam_index = donor_SortSampleBam.output_bam_index,
        output_bam_basename = donor_recalibrated_bam_basename,
        recalibration_report = donor_GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        # We need disk to localize the sharded bam and the sharded output due to the scatter.
        disk_size = ceil((donor_agg_bam_size * 3) + ref_size) + disk_pad,
        compression_level = compression_level,
        preemptible_tries = agg_preemptible_tries
    }
  }

  call GatherBamFiles as donor_GatherBamFiles{
    input:
      input_bams = donor_ApplyBQSR.recalibrated_bam,
      output_bam_basename = donor_output_basename,
      disk_size = (2 * donor_agg_bam_size) + disk_pad,
      compression_level = compression_level,
      preemptible_tries = agg_preemptible_tries
  }

  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call BaseRecalibrator as host_BaseRecalibrator {
      input:
        input_bam = host_SortSampleBam.output_bam,
        input_bam_index = host_SortSampleBam.output_bam_index,
        recalibration_report_filename = host_output_basename + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        disk_size = ceil(host_agg_bam_size + ref_size + dbsnp_size) + disk_pad,
        preemptible_tries = agg_preemptible_tries
    }
  }

  # Merge the recalibration reports resulting from by-interval recalibration
  # The reports are always the same size
  call GatherBqsrReports as host_GatherBqsrReports {
    input:
      input_bqsr_reports = host_BaseRecalibrator.recalibration_report,
      output_report_filename = host_output_basename + ".recal_data.csv",
      disk_size = disk_pad,
      preemptible_tries = preemptible_tries
  }

  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {
    # Apply the recalibration model by interval
    call ApplyBQSR as host_ApplyBQSR {
      input:
        input_bam = host_SortSampleBam.output_bam,
        input_bam_index = host_SortSampleBam.output_bam_index,
        output_bam_basename = host_recalibrated_bam_basename,
        recalibration_report = host_GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        # We need disk to localize the sharded bam and the sharded output due to the scatter.
        disk_size = ceil((host_agg_bam_size * 3) + ref_size) + disk_pad,
        compression_level = compression_level,
        preemptible_tries = agg_preemptible_tries
    }
  }

  call GatherBamFiles as host_GatherBamFiles {
    input:
      input_bams = host_ApplyBQSR.recalibrated_bam,
      output_bam_basename = donor_output_basename,
      disk_size = (2 * host_agg_bam_size) + disk_pad,
      compression_level = compression_level,
      preemptible_tries = agg_preemptible_tries
  }

  call interval_list_to_bed {
      input:
          interval_list=capture_interval_list
  }

  call deep_variant as donor_deepvariant {
    input: 
      sample=donor_output_basename,
      capture_bed=interval_list_to_bed.bed,
      bam=donor_GatherBamFiles.output_bam,
      bai=donor_GatherBamFiles.output_bam_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      runtime_docker = runtime_docker
  }

    call bgzip as donor_bgzip {
        input:
            sample=donor_output_basename,
            uncompressed_vcf=donor_deepvariant.vcf,
            runtime_disk = runtime_disk
    }

  call deep_variant as host_deepvariant {
    input: 
      sample=host_output_basename,
      capture_bed=interval_list_to_bed.bed,
      bam=host_GatherBamFiles.output_bam,
      bai=host_GatherBamFiles.output_bam_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      runtime_docker = runtime_docker
  }

    call bgzip as host_bgzip {
        input:
            sample=host_output_basename,
            uncompressed_vcf=host_deepvariant.vcf,
            runtime_disk = runtime_disk
    }

    call bcf_filter {
            input:
              donor_input_vcf = donor_bgzip.filtered_vcf,
              host_input_vcf = host_bgzip.filtered_vcf
    }

    # Donor annotating 
    call Funcotate as FuncotateDonor {
            input:
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    ref_dict = ref_dict,
                    input_vcf = bcf_filter.donor_vcf,
                    input_vcf_idx = bcf_filter.donor_vcf_tbi,
                    sample = donor_output_basename,

                    data_sources_tar_gz = data_sources_tar_gz
    }

    # Host annotating 
    call Funcotate as FuncotateHost {
            input:
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    ref_dict = ref_dict,
                    input_vcf = bcf_filter.host_vcf,
                    input_vcf_idx = bcf_filter.host_vcf_tbi,
                    sample = host_output_basename,

                    data_sources_tar_gz = data_sources_tar_gz
    }

    call bmt {
            input: 
                    donor_input_vcf = FuncotateDonor.funcotated_output_file,
                    host_input_vcf = FuncotateHost.funcotated_output_file,
                    tissue = tissue,
                    donor_sex = donor_sex,
                    host_sex = host_sex,
                    gvhd_filter = gvhd_filter
    }

    #Array[File] discordances_file = select_first([bmt.specific_discordances, bmt.gvhd_specific_discordances])

    # Donor specific_discordances = tissue_specific_discordances.txt
    call blastFileGenerationForBlast as gvl_DonorblastFileGenerationForBlast {
            input:
                    specific_discordances = bmt.specific_discordances
    }

    # Creating donor custom proteome using makeblastdb 
    # donor custom proteome against donor peptides as query 
    call blastOperations as gvl_blastOperationsDonor {
            input: 
                    proteome_filename = bmt.donor_custom_proteome,
                    makeblastdb_out = "donor_db",
                    peptide_file = gvl_DonorblastFileGenerationForBlast.peptides_for_blastp, # for blastp
                    blastp_out = donor_output_basename+"_blastp_peptides_out.csv",
                    sample = donor_output_basename
    }

    # For donor : specific_discordances = tissue_specific_discordances.txt -> donor_discordances_after_blast.txt
    call blastOperationsPostProcessing as gvl_DonorBlastOperationsPostProcessing {
            input: 
                    specific_discordances = bmt.specific_discordances,
                    input_blastp_result = gvl_blastOperationsDonor.blastp_result,
                    sample = donor_output_basename,
                    sample_type = "donor",
                    condition="GvL"
    }


    # Donor discordances after blast as input 
    call blastFileGenerationForBlast as gvl_HostblastFileGenerationForBlast {
            input:
                    specific_discordances = gvl_DonorBlastOperationsPostProcessing.discordances_after_blast
    }

    # Creating host custom proteome using makeblastdb 
    # host custom proteome against host peptides as query 
    call blastOperations as gvl_blastOperationsHost {
            input: 
                    proteome_filename = bmt.host_custom_proteome,
                    makeblastdb_out = "host_db",
                    peptide_file = gvl_HostblastFileGenerationForBlast.peptides_for_blastp, # for blastp 
                    blastp_out = host_output_basename+"_blastp_peptides_out.csv",
                    sample = host_output_basename
    }

    # For host : specific_discordances = donor_discordances_after_blast.txt with genes to remove -> host_discordances_after_blast.txt
    call blastOperationsPostProcessing as gvl_HostBlastOperationsPostProcessing {
            input: 
                    specific_discordances = gvl_DonorBlastOperationsPostProcessing.discordances_after_blast,
                    input_blastp_result = gvl_blastOperationsHost.blastp_result,
                    sample_type = "host",
                    sample = host_output_basename,
                    condition="GvL"
    }

    call HLAthena_preprocessing as gvl_HLAthena_preprocessing {
            input:
                    discordances_after_blast = gvl_HostBlastOperationsPostProcessing.allgenes_pumas
    }

    call HLAthena.parse_alleles_task {
            input : 
                    preemptible = 3,
                    memoryGB = 8,
                    diskGB = 10,

                    alleles_file = hlatypes_file
    }


    Array[String] gvl_alleles_parsed = read_lines(parse_alleles_task.alleles_parsed)
    if (aggregate_pep) {
            call HLAthena.aggregate_pep_task as gvl_aggregate_pep_task {
                    input:
                            preemptible = 3,
                            memoryGB = 8,
                            diskGB = 10,

                            peptide_list = gvl_HLAthena_preprocessing.HLAthena_preprocessed_peptides,
                            peptide_col_name = peptide_col_name,
                            exists_expr = exists_expr,
                            expr_col_name = expr_col_name
            }
    }

    ### Split merged peptides by length
    call HLAthena.split_peptides_len_task as gvl_split_peptides_len_task {
            input:
                preemptible = 3,
                memoryGB = 8,
                diskGB = 10,

                peptide_list = gvl_aggregate_pep_task.peptide_list_agg,
                peptide_col_name = peptide_col_name,
                lens = lens
    }

    Array[String] gvl_lens_present = read_lines(gvl_split_peptides_len_task.lens_present)

    ### Generate peptide sequences features: dummy encoding only 
    ### (blosum and fuzzy generated on demand from the dummy encoding upon prediction)
    scatter (pepfilelen in gvl_split_peptides_len_task.pep_len_files) {
            call HLAthena.featurize_encoding_task as gvl_featurize_encoding_task {
                    input:
                            preemptible = 3,
                            memoryGB = 8,
                            diskGB = 10,

                            peptide_len_list = pepfilelen,
                            encoding = "dummy",
                            peptide_col_name = peptide_col_name
            }
    }

    ### Predict with MS models (now includes ranks)
    scatter (allele_featfile in cross(gvl_alleles_parsed, gvl_featurize_encoding_task.feature_files)) {
            call HLAthena.predict_ms_allele_task as gvl_predict_ms_allele_task {
                    input:
                            preemptible = 3,
                            memoryGB = 8,
                            diskGB = 10,

                            patient = sample_id,
                            featfile = allele_featfile.right,
                            lens = gvl_lens_present,
                            exists_ctex = exists_ctex,
                            exists_expr = exists_expr,
                            allele = allele_featfile.left,
                            models_path = "gs://msmodels/",
                            peptide_col_name = peptide_col_name,
                            features_all = "features_AAPos_AAPCA_LogTPM_CNN_Kidera_Gene",
                            feature_sets = "features_AAPos_AAPCA_Kidera"
            }
    }

    ### Merge predictions into ensemble model scores: MSIntrinsic, MSIntrinsicC, MSIntrinsicEC; optinally merge NetMHC scores
    scatter (len in gvl_lens_present) {
            call HLAthena.merge_mspreds_task as gvl_merge_mspreds_task {
                    input:
                            preemptible = 3,
                            memoryGB_merge = 16,
                            #diskGB_merge = 20, - PRESENT IN THE TASK BUT ERROR THAT IT IS NOT PRESENT 
                            diskGB = 10,
                            
                            patient = sample_id,
                            alleles = gvl_alleles_parsed,
                            len = len,
                            exists_ctex = exists_ctex,
                            exists_expr = exists_expr,
                            mspred_files = gvl_predict_ms_allele_task.mspred_files,
                            peptide_col_name = peptide_col_name
                            #run_netmhc = run_netmhc,
                            #netmhc_EL = predict_netmhc_task.netmhc_EL,
                            #netmhc_BA = predict_netmhc_task.netmhc_BA
            }
    }

    ### Convert scores to ranks
    scatter (preds_file in gvl_merge_mspreds_task.mspreds_file_wide) {
            call HLAthena.get_ranks_task as gvl_get_ranks_task {
                    input:
                            preemptible = 3,
                            memoryGB_ranks = 16,
                            diskGB = 10,

                            patient = sample_id,
                            alleles = gvl_alleles_parsed,
                            lens = gvl_lens_present,
                            exists_ctex = exists_ctex,
                            exists_expr = exists_expr,
                            #run_netmhc = run_netmhc,
                            models_path = "gs://msmodels/",
                            peptide_col_name = peptide_col_name,
                            preds_file = preds_file
            }
    }

    ### Concatenate ranks files for all lengths
    call HLAthena.ranks_concat_lens_task as gvl_ranks_concat_lens_task {
            input:
                    preemptible = 3,
                    memoryGB_assign = 16,
                    #diskGB_ranks_concat = 20,
                    diskGB = 10,

                    patient = sample_id,
                    peptide_col_name = peptide_col_name,
                    exists_ctex = exists_ctex,
                    exists_expr = exists_expr,
                    #run_netmhc = run_netmhc,
                    assign_by_ranks_or_scores = "ranks",
                    assign_threshold = "0.5",
                    assign_colors = "#8c0a82, #ab0c9f, #048edb, #00a0fa, #5F36A0, #ff6000, darkorange, #e05400, grey",
                    ranks_len_files = gvl_get_ranks_task.ranks_file
    }

    call HLAthena_postprocessing as gvl_HLAthena_postprocessing {
            input :
                    sample_predictions = gvl_ranks_concat_lens_task.sample_predictions,
                    discordances_after_blast = gvl_HostBlastOperationsPostProcessing.allgenes_pumas,
                    tissue = tissue,
                    sample_name = sample_id + "_GvL"
    }



    ############ GvHD ############


    # Donor specific_discordances = tissue_specific_discordances.txt
    call blastFileGenerationForBlast as gvhd_DonorblastFileGenerationForBlast {
            input:
                    specific_discordances = bmt.gvhd_specific_discordances
    }

    # Creating donor custom proteome using makeblastdb 
    # donor custom proteome against donor peptides as query 
    call blastOperations as gvhd_blastOperationsDonor {
            input: 
                    proteome_filename = bmt.donor_custom_proteome,
                    makeblastdb_out = "donor_db",
                    peptide_file = gvhd_DonorblastFileGenerationForBlast.peptides_for_blastp, # for blastp
                    blastp_out = donor_output_basename+"_blastp_peptides_out.csv",
                    sample = donor_output_basename
    }

    # For donor : specific_discordances = tissue_specific_discordances.txt -> donor_discordances_after_blast.txt
    call blastOperationsPostProcessing as gvhd_DonorBlastOperationsPostProcessing {
            input: 
                    specific_discordances = bmt.gvhd_specific_discordances,
                    input_blastp_result = gvhd_blastOperationsDonor.blastp_result,
                    sample = donor_output_basename,
                    sample_type = "donor",
                    condition="GvHD"
    }


    # Donor discordances after blast as input 
    call blastFileGenerationForBlast as gvhd_HostblastFileGenerationForBlast {
            input:
                    specific_discordances = gvhd_DonorBlastOperationsPostProcessing.discordances_after_blast
    }

    # Creating host custom proteome using makeblastdb 
    # host custom proteome against host peptides as query 
    call blastOperations as gvhd_blastOperationsHost {
            input: 
                    proteome_filename = bmt.host_custom_proteome,
                    makeblastdb_out = "host_db",
                    peptide_file = gvhd_HostblastFileGenerationForBlast.peptides_for_blastp, # for blastp 
                    blastp_out = host_output_basename+"_blastp_peptides_out.csv",
                    sample = host_output_basename
    }

    # For host : specific_discordances = donor_discordances_after_blast.txt with genes to remove -> host_discordances_after_blast.txt
    call blastOperationsPostProcessing as gvhd_HostBlastOperationsPostProcessing {
            input: 
                    specific_discordances = gvhd_DonorBlastOperationsPostProcessing.discordances_after_blast,
                    input_blastp_result = gvhd_blastOperationsHost.blastp_result,
                    sample_type = "host",
                    sample = host_output_basename,
                    condition="GvHD"
    }

    call HLAthena_preprocessing as gvhd_HLAthena_preprocessing {
            input:
                    discordances_after_blast = gvhd_HostBlastOperationsPostProcessing.allgenes_pumas
    }

    Array[String] gvhd_alleles_parsed = read_lines(parse_alleles_task.alleles_parsed)
    if (aggregate_pep) {
            call HLAthena.aggregate_pep_task as gvhd_aggregate_pep_task {
                    input:
                            preemptible = 3,
                            memoryGB = 8,
                            diskGB = 10,

                            peptide_list = gvhd_HLAthena_preprocessing.HLAthena_preprocessed_peptides,
                            peptide_col_name = peptide_col_name,
                            exists_expr = exists_expr,
                            expr_col_name = expr_col_name
            }
    }

    ### Split merged peptides by length
    call HLAthena.split_peptides_len_task as gvhd_split_peptides_len_task {
            input:
                preemptible = 3,
                memoryGB = 8,
                diskGB = 10,

                peptide_list = gvhd_aggregate_pep_task.peptide_list_agg,
                peptide_col_name = peptide_col_name,
                lens = lens
    }

    Array[String] gvhd_lens_present = read_lines(gvhd_split_peptides_len_task.lens_present)

    ### Generate peptide sequences features: dummy encoding only 
    ### (blosum and fuzzy generated on demand from the dummy encoding upon prediction)
    scatter (pepfilelen in gvhd_split_peptides_len_task.pep_len_files) {
            call HLAthena.featurize_encoding_task as gvhd_featurize_encoding_task {
                    input:
                            preemptible = 3,
                            memoryGB = 8,
                            diskGB = 10,

                            peptide_len_list = pepfilelen,
                            encoding = "dummy",
                            peptide_col_name = peptide_col_name
            }
    }

    ### Predict with MS models (now includes ranks)
    scatter (allele_featfile in cross(gvhd_alleles_parsed, gvhd_featurize_encoding_task.feature_files)) {
            call HLAthena.predict_ms_allele_task as gvhd_predict_ms_allele_task {
                    input:
                            preemptible = 3,
                            memoryGB = 8,
                            diskGB = 10,

                            patient = sample_id,
                            featfile = allele_featfile.right,
                            lens = gvhd_lens_present,
                            exists_ctex = exists_ctex,
                            exists_expr = exists_expr,
                            allele = allele_featfile.left,
                            models_path = "gs://msmodels/",
                            peptide_col_name = peptide_col_name,
                            features_all = "features_AAPos_AAPCA_LogTPM_CNN_Kidera_Gene",
                            feature_sets = "features_AAPos_AAPCA_Kidera"
            }
    }

    ### Merge predictions into ensemble model scores: MSIntrinsic, MSIntrinsicC, MSIntrinsicEC; optinally merge NetMHC scores
    scatter (len in gvhd_lens_present) {
            call HLAthena.merge_mspreds_task as gvhd_merge_mspreds_task {
                    input:
                            preemptible = 3,
                            memoryGB_merge = 16,
                            #diskGB_merge = 20, - PRESENT IN THE TASK BUT ERROR THAT IT IS NOT PRESENT 
                            diskGB = 10,
                            
                            patient = sample_id,
                            alleles = gvhd_alleles_parsed,
                            len = len,
                            exists_ctex = exists_ctex,
                            exists_expr = exists_expr,
                            mspred_files = gvhd_predict_ms_allele_task.mspred_files,
                            peptide_col_name = peptide_col_name
                            #run_netmhc = run_netmhc,
                            #netmhc_EL = predict_netmhc_task.netmhc_EL,
                            #netmhc_BA = predict_netmhc_task.netmhc_BA
            }
    }

    ### Convert scores to ranks
    scatter (preds_file in gvhd_merge_mspreds_task.mspreds_file_wide) {
            call HLAthena.get_ranks_task as gvhd_get_ranks_task {
                    input:
                            preemptible = 3,
                            memoryGB_ranks = 16,
                            diskGB = 10,

                            patient = sample_id,
                            alleles = gvhd_alleles_parsed,
                            lens = gvhd_lens_present,
                            exists_ctex = exists_ctex,
                            exists_expr = exists_expr,
                            #run_netmhc = run_netmhc,
                            models_path = "gs://msmodels/",
                            peptide_col_name = peptide_col_name,
                            preds_file = preds_file
            }
    }

    ### Concatenate ranks files for all lengths
    call HLAthena.ranks_concat_lens_task as gvhd_ranks_concat_lens_task {
            input:
                    preemptible = 3,
                    memoryGB_assign = 16,
                    #diskGB_ranks_concat = 20,
                    diskGB = 10,

                    patient = sample_id,
                    peptide_col_name = peptide_col_name,
                    exists_ctex = exists_ctex,
                    exists_expr = exists_expr,
                    #run_netmhc = run_netmhc,
                    assign_by_ranks_or_scores = "ranks",
                    assign_threshold = "0.5",
                    assign_colors = "#8c0a82, #ab0c9f, #048edb, #00a0fa, #5F36A0, #ff6000, darkorange, #e05400, grey",
                    ranks_len_files = gvhd_get_ranks_task.ranks_file
    }

    call HLAthena_postprocessing as gvhd_HLAthena_postprocessing {
            input :
                    sample_predictions = gvhd_ranks_concat_lens_task.sample_predictions,
                    discordances_after_blast = gvhd_HostBlastOperationsPostProcessing.allgenes_pumas,
                    tissue = tissue,
                    sample_name = sample_id + "_GvHD"
    }

    call merge_summary_csv_files {
      input:
        bcf_dict = bcf_filter.bcf_summary,
        bmt_gvl_dict = bmt.GvL_dict_nums,
        gvl_postBlast_dict = gvl_HostBlastOperationsPostProcessing.GvL_postBlast_dict_nums,
        gvl_postHLAthena_dict = gvl_HLAthena_postprocessing.GvL_postHLAthena_dict_nums,
        bmt_gvhd_dict = bmt.GvHD_dict_nums,
        gvhd_postBlast_dict = gvhd_HostBlastOperationsPostProcessing.GvHD_postBlast_dict_nums,
        gvhd_postHLAthena_dict = gvhd_HLAthena_postprocessing.GvHD_postHLAthena_dict_nums

    }

    output {
      # BMT OUTPUT 
      File specific_discordances_file = bmt.specific_discordances 
      File gvhd_specific_discordances_file = bmt.gvhd_specific_discordances 

      # HLAthena OUTPUT
      File gvl_sample_predictions_file = gvl_ranks_concat_lens_task.sample_predictions
      File gvhd_sample_predictions_file = gvhd_ranks_concat_lens_task.sample_predictions

      # HLAthena POST PROCESSING OUTPUT
      File gvl_binding_putativeMinorAntigens_file = gvl_HLAthena_postprocessing.binding_putativeMinorAntigens
      File gvhd_binding_putativeMinorAntigens_file = gvhd_HLAthena_postprocessing.binding_putativeMinorAntigens

      # FINAL OUTPUT
      File pipeline_run_summary_csv = merge_summary_csv_files.output_summary_file_result
      File pipeline_run_summary_csv_2 = merge_summary_csv_files.output_summary_file_result_2
    }
}


# REALIGNING FROM HG38 TO HG19 

# Conversion from CRAM to BAM with intermediate SAM file of hg38 genome built 
task CramToBam {
  File ref_fasta
  File ref_fasta_index
  File cram_file
  String output_basename

  Int disk_size

command <<<
  set -e
  set -o pipefail

  samtools view -h -T ${ref_fasta} ${cram_file} |
  samtools view -b -o ${output_basename}.bam -
  samtools index -b ${output_basename}.bam
  mv ${output_basename}.bam.bai ${output_basename}.bai
  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "${output_basename}.bam"
    File output_bam_index = "${output_basename}.bai"
  }
}



# Conversion of the BAM into uBAM
task RevertSam {
  File input_bam
  Int disk_size
  String output_basename

  command {
    java -Xmx1000m -jar /usr/gitc/picard.jar \
    RevertSam \
    INPUT=${input_bam} \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=LENIENT \
    ATTRIBUTE_TO_CLEAR=FT \
    ATTRIBUTE_TO_CLEAR=CO \
    OUTPUT=${output_basename}.bam
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
    disks: "local-disk " + disk_size + " HDD"
    memory: "1200 MB"
    cpus: 8
  }
  output {
    File unmapped_bam = "${output_basename}.bam"
  }
}

# Sorts the uBAM wrt "queryname"
# queryname: every read has a queryname, first of the description. This sorts according to the read groups, but doesn't necessarily sort the alignments (reads)
# queryname sorting is important in order to run MergeBamAlignment in the next task (GATK recommended)

task SortSam {
  File input_bam
  String sorted_bam_name
  Int disk_size

  command {
    java -Xmx3000m -jar /usr/gitc/picard.jar \
    SortSam \
    INPUT=${input_bam} \
    OUTPUT=${sorted_bam_name} \
    SORT_ORDER=queryname \
    MAX_RECORDS_IN_RAM=1000000 # to avoid spillage onto disk in case of low memory (sorting happens in the memory and is written onto the disk)
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
    disks: "local-disk " + disk_size + " HDD"
    bootDiskSizeGb: 12
    memory: "20 GB"
    preemptible: 3
  }
  output {
    File sorted_bam = "${sorted_bam_name}"
  }
}



# 3 tasks 
# 1. Sam/Bam to 2 Fastq (paired-end)
# 2. 2 Fastq, hg19 ref file and bwa index files for BWA aligner to align -> BAM file
# 3. aligned BAM file is merged with uBAM (to include metadata: sample alias, library, barcodes, etc) which is not present in fastq files. Yields an unsorted BAM.

task SamToFastqAndBwaMemAndMba {
  File input_bam
  String bwa_commandline
  String bwa_version
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
  # listing the reference contigs that are "alternative".

  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
  Float disk_size
  Int compression_level
  Int preemptible_tries

  command <<<
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
      java -Xms5000m -jar /usr/gitc/picard.jar \
        SamToFastq \
        INPUT=${input_bam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true | \
      /usr/gitc/${bwa_commandline} /dev/stdin - 2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) | \
      java -Dsamjdk.compression_level=${compression_level} -Xms3000m -jar /usr/gitc/picard.jar \
        MergeBamAlignment \
        VALIDATION_STRINGENCY=SILENT \
        EXPECTED_ORIENTATIONS=FR \
        ATTRIBUTES_TO_RETAIN=X0 \
        ATTRIBUTES_TO_REMOVE=NM \
        ATTRIBUTES_TO_REMOVE=MD \
        ALIGNED_BAM=/dev/stdin \
        UNMAPPED_BAM=${input_bam} \
        OUTPUT=${output_bam_basename}.bam \
        REFERENCE_SEQUENCE=${ref_fasta} \
        PAIRED_RUN=true \
        SORT_ORDER="unsorted" \
        IS_BISULFITE_SEQUENCE=false \
        ALIGNED_READS_ONLY=false \
        CLIP_ADAPTERS=false \
        MAX_RECORDS_IN_RAM=2000000 \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        PROGRAM_RECORD_ID="bwamem" \
        PROGRAM_GROUP_VERSION="${bwa_version}" \
        PROGRAM_GROUP_COMMAND_LINE="${bwa_commandline}" \
        PROGRAM_GROUP_NAME="bwamem" \
        ALIGNER_PROPER_PAIR_FLAGS=true \
        UNMAP_CONTAMINANT_READS=true \
        ADD_PG_TAG_TO_READS=false

  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "14 GB"
    bootDiskSizeGb: 12
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
    cpu: "16"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
  }
}



# Marks duplicate reads (defined as originating from a single fragment of DNA)
# Compares sequences in the 5 prime positions of both reads 
# Differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums of their base-quality scores

task MarkDuplicates {
  File input_bams
  String output_bam_basename
  String metrics_filename
  Float disk_size
  Int compression_level
  Int preemptible_tries

  # The program default for READ_NAME_REGEX is appropriate in nearly every case.
  # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
  # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
  String? read_name_regex

 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
 # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
    java -Dsamjdk.compression_level=${compression_level} -Xms4000m -jar /usr/gitc/picard.jar \
      MarkDuplicates \
      INPUT=${input_bams} \
      OUTPUT=${output_bam_basename}.bam \
      METRICS_FILE=${metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      ${"READ_NAME_REGEX=" + read_name_regex} \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname" \
      ADD_PG_TAG_TO_READS=false
  }
  runtime {
    preemptible: preemptible_tries
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
    memory: "7 GB"
    bootDiskSizeGb: 12
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File duplicate_metrics = "${metrics_filename}"
  }
}


# Sorts wrt "coordinate" 
# 2 step process: 1. sorts the RNAME in @SQ/SN (header), 2. sorts subgroups based on POS (left-most read is mapped and so on)
task SortSampleBam {
  File input_bam
  String output_bam_basename
  Int preemptible_tries
  Int compression_level
  Float disk_size

  command {
    java -Dsamjdk.compression_level=${compression_level} -Xms4000m -jar /usr/gitc/picard.jar \
      SortSam \
      INPUT=${input_bam} \
      OUTPUT=${output_bam_basename}.bam \
      SORT_ORDER="coordinate" \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true \
      MAX_RECORDS_IN_RAM=300000

  }
  runtime {
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
    cpu: "1"
    bootDiskSizeGb: 12
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
    memory: "20 GB"
    preemptible: preemptible_tries
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}



# Creating genomic intervals over chromosomes 
# Split based on the length of the sequence in reference sequence - The longest sequence is the threshold, sequence lengths madding up to longest sequences are binned together 
task CreateSequenceGroupingTSV {
  File ref_dict
  Int preemptible_tries

  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter.
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python <<CODE
    with open("${ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
    # the last element after a :, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>
  runtime {
    memory: "2 GB"
    docker: "python:2.7"
    preemptible: preemptible_tries
  }
  output {
    Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}



# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  File input_bam
  File input_bam_index
  String recalibration_report_filename
  Array[String] sequence_group_interval
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Float disk_size
  Int preemptible_tries

  command {
    /usr/gitc/gatk4/gatk-launch --javaOptions "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log -Xms4000m" \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --useOriginalQualities \
      -O ${recalibration_report_filename} \
      -knownSites ${dbSNP_vcf} \
      -knownSites ${sep=" -knownSites " known_indels_sites_VCFs} \
      -L ${sep=" -L " sequence_group_interval}
  }
  runtime {
    preemptible: preemptible_tries
    memory: "6 GB"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
  }
  output {
    File recalibration_report = "${recalibration_report_filename}"
  }
}


# Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
  Array[File] input_bqsr_reports
  String output_report_filename
  Int disk_size
  Int preemptible_tries

  command {
    /usr/gitc/gatk4/gatk-launch --javaOptions "-Xms9000m" \
      GatherBQSRReports \
      -I ${sep=' -I ' input_bqsr_reports} \
      -O ${output_report_filename}
    }
  runtime {
    preemptible: preemptible_tries
    memory: "10 GB"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bqsr_report = "${output_report_filename}"
  }
}


# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  File input_bam
  File input_bam_index
  String output_bam_basename
  File recalibration_report
  Array[String] sequence_group_interval
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Float disk_size
  Int compression_level
  Int preemptible_tries

  command {
   /usr/gitc/gatk4/gatk-launch --javaOptions "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
      -XX:+PrintGCDetails -Xloggc:gc_log.log \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level=${compression_level} -Xms9000m" \
      ApplyBQSR \
      --createOutputBamMD5 \
      --addOutputSAMProgramRecord \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --useOriginalQualities \
      -O ${output_bam_basename}.bam \
      -bqsr ${recalibration_report} \
      -L ${sep=" -L " sequence_group_interval}
  }
  runtime {
    preemptible: preemptible_tries
    memory: "10 GB"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
  }
  output {
    File recalibrated_bam = "${output_bam_basename}.bam"
    File recalibrated_bam_checksum = "${output_bam_basename}.bam.md5"
  }
}


# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
  Array[File] input_bams
  String output_bam_basename
  Float disk_size
  Int compression_level
  Int preemptible_tries

  command {
    java -Dsamjdk.compression_level=${compression_level} -Xms2000m -jar /usr/gitc/picard.jar \
      GatherBamFiles \
      INPUT=${sep=' INPUT=' input_bams} \
      OUTPUT=${output_bam_basename}.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true
    }
  runtime {
    preemptible: preemptible_tries
    memory: "3 GB"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}


# Interval list (coding regions of hg19) is converted into bed format 
task interval_list_to_bed {
    File interval_list
    String bed_path = sub(basename(interval_list), 'interval_list', 'bed')

    command <<<
    set -xeuo pipefail

    # interval lists have headers that need to be removed and are 1-indexed
    # see also https://www.biostars.org/p/84686/
    grep -v '^@' ${interval_list} \
    | awk -v OFS='\t' '{print $1, $2 - 1, $3}' \
    | sort -k1,1 -k2,2n -k3,3n \
    | bedtools merge \
    > ${bed_path}
    >>>

    output {
        File bed = '${bed_path}'
    }

    runtime {
        memory: '1 GB'
        disks: 'local-disk 1 HDD'
        preemptible: 3
        docker: 'quay.io/biocontainers/bedtools:2.28.0--hdf88d34_0'
    }
}



# Germline Variant calling 
task deep_variant {
    String sample
    File bam
    File bai
    String model_type = 'WES'
    File ref_fasta
    File ref_fasta_index
    File? capture_bed

    Int runtime_cpus = 32 
    String runtime_docker
    Int runtime_memory = ceil(1.1 * runtime_cpus)
    Int runtime_disk_buffer = 3
    Int runtime_disk = ceil(1.15 * (size(ref_fasta, 'G') + size(bam, 'G')) + runtime_disk_buffer)
    Int runtime_preemptible = 3
    Int resource_log_interval = 10 

    command <<<
    # log resource usage for debugging purposes
    function runtimeInfo() {
        echo [$(date)]
        echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
        echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
        echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
    }
    while true;
        do runtimeInfo >> resource.log;
        sleep ${resource_log_interval};
    done &
    lscpu

    set -xeuo pipefail

    # make symbolic links to ensure BAM and index are in expected structure even after localization
    ln -s ${bai} reads.bai
    ln -s ${bam} reads.bam

    # make symbolic links to ensure reference and index are in expected structure even after localization
    ln -s ${ref_fasta} reference.fa
    ln -s ${ref_fasta_index} reference.fa.fai

    mkdir deepvariant_tmp

    /opt/deepvariant/bin/run_deepvariant \
        --model_type=${model_type} \
        --ref=reference.fa \
        --reads=reads.bam \
        --regions=${capture_bed} \
        --intermediate_results_dir=deepvariant_tmp \
        --output_vcf=${sample}.vcf \
        --num_shards=${runtime_cpus}
    >>>

    output {
        File vcf = '${sample}.vcf'
        File resource_log = 'resource.log'
    }
    
    runtime {
        memory: '${runtime_memory} GB'
        cpu: '${runtime_cpus}'
        disks: 'local-disk ${runtime_disk} SSD'
        preemptible: '${runtime_preemptible}'
        docker: '${runtime_docker}'
    }
}


# gunzips the VCF files from deepvariant and contains only "PASS" variants 
task bgzip {
    String sample
    File uncompressed_vcf
    Int runtime_disk

    command <<<
    set -xeuo pipefail

    bcftools view -Oz -o ${sample}.vcf.gz ${uncompressed_vcf}
    bcftools index --tbi ${sample}.vcf.gz

    # create version of VCF with only PASSing variants
    bcftools view -Oz -o ${sample}_filtered.vcf.gz -f PASS ${uncompressed_vcf}
    bcftools index --tbi ${sample}_filtered.vcf.gz
    >>>

    output {
        File vcf = '${sample}.vcf.gz'
        File vcf_index = '${sample}.vcf.gz.tbi'
        File filtered_vcf = '${sample}_filtered.vcf.gz'
        File filtered_vcf_index = '${sample}_filtered.vcf.gz.tbi'
    }

    runtime {
        memory: '1 GB'
        disks: 'local-disk ${runtime_disk} HDD'
        preemptible: 3
        docker: 'quay.io/biocontainers/bcftools:1.9--ha228f0b_3'
    }
}


# Filters "PASS" and "MIS" variants 
task bcf_filter {

  # Filtering the Deepvariant VCF for donor and host files to remore ref_calls and mis 

        File donor_input_vcf
        File host_input_vcf

        command {
                bcftools filter ${donor_input_vcf} -i 'FILTER="PASS" && GT!="mis"' > donor_bcf_filtered.vcf.gz -O z 
                tabix -p vcf donor_bcf_filtered.vcf.gz

                bcftools filter ${host_input_vcf} -i 'FILTER="PASS" && GT!="mis"' > host_bcf_filtered.vcf.gz -O z
                tabix -p vcf host_bcf_filtered.vcf.gz

                donor_bcf_num=$(bcftools filter ${donor_input_vcf} -i 'FILTER="PASS" && GT!="mis"' | grep -c '^[^#]')
                host_bcf_num=$(bcftools filter ${host_input_vcf} -i 'FILTER="PASS" && GT!="mis"' | grep -c '^[^#]')

                echo "host_VCF_rows","donor_VCF_rows" > bcf_summary.csv
                echo "$host_bcf_num","$donor_bcf_num" >> bcf_summary.csv

                rm -r ${donor_input_vcf}
                rm -r ${host_input_vcf}             
        }

        output {
                File donor_vcf="donor_bcf_filtered.vcf.gz"
                File donor_vcf_tbi="donor_bcf_filtered.vcf.gz.tbi"

                File host_vcf="host_bcf_filtered.vcf.gz"
                File host_vcf_tbi="host_bcf_filtered.vcf.gz.tbi"

                File bcf_summary="bcf_summary.csv"
        }

        runtime {
          docker : "dceoy/bcftools"
        }

        meta {
                author: "Nidhi Hookeri"
                email: "nhookeri@broadinstitute.org"
                description: "Filters the vcf file."
        }
}


# Functional Annotator 
task Funcotate {
    
  # Annotating the BCF filtered VCF files for donor and host

        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File input_vcf
        File input_vcf_idx
        #String reference_version
        #Boolean compress
        String sample

        #Int num_threads

        File? data_sources_tar_gz = "gs://broad-public-datasets/funcotator/funcotator_dataSources.v1.6.20190124g.tar.gz"
        #String? transcript_selection_mode = "ALL"

        String output_vcf = sample + "_funcotated" + ".vcf"
        String output_vcf_idx = output_vcf + ".idx"
        
        command <<<
                # Extract the tar.gz:
                echo "Extracting data sources tar/gzip file..."
                mkdir datasources_dir
                tar zxvf ${data_sources_tar_gz} -C datasources_dir --strip-components 1
                DATA_SOURCES_FOLDER="$PWD/datasources_dir"

                # Annotating vcf file
                gatk --java-options "-Xmx2048m" Funcotator \
                     --data-sources-path $DATA_SOURCES_FOLDER \
                     --ref-version "hg19" \
                     --output-file-format "VCF" \
                     -R ${ref_fasta} \
                     -V ${input_vcf} \
                     -O ${output_vcf} \
                     --transcript-selection-mode "ALL" \
                     --create-output-variant-index
                bgzip ${output_vcf}

                rm -r ${input_vcf}
                rm -r ${input_vcf_idx}
                rm -r ${ref_fasta}
                rm -r ${ref_fasta_index}
                rm -r ${ref_dict}
                rm -r ${data_sources_tar_gz}
        >>>

        runtime {
                docker: "us.gcr.io/broad-gatk/gatk:4.1.4.0"
                cpu: 32
        }

        output {
                File funcotated_output_file = "${output_vcf}.gz"
                File funcotated_output_file_index = "${output_vcf_idx}"
        }
}



task bmt {
        File donor_input_vcf
        File host_input_vcf
        String tissue
        String donor_sex
        String host_sex
        String gvhd_filter

        command {
                pwd
                df
                python3 /pipeline/bmt-simulation.py -donor_vcf ${donor_input_vcf} -host_vcf ${host_input_vcf} -donor_sex ${donor_sex} -host_sex ${host_sex} -tissue ${tissue} -gvhd ${gvhd_filter}
        }

        runtime {
                docker : "nidhihookeri/minors-pipeline-2"
                cpu: 32
        }
        output {
                String result = stdout()
                #File donor_annotation = "donor_annotations.txt"
                #File donor_not_mutated = "donor_not_mutated.txt"

                #File host_annotation = "host_annotations.txt"
                #File host_not_mutated = "host_not_mutated.txt"

                #File donor_annotations_withprot = "donor_annotations_withprot.txt"
                #File donor_mutated_proteins = "donor_mutated_proteins.txt"
                #File donor_mutation_indices = "donor_mutation_indices.txt"

                #File host_annotations_withprot = "host_annotations_withprot.txt"
                #File host_mutated_proteins = "host_mutated_proteins.txt"
                #File host_mutation_indices = "host_mutation_indices.txt"

                #File donor_and_host_data_merged = "donor_and_host_data_merged.txt"

                File donor_custom_proteome = "donor_custom_proteome.fasta"
                File host_custom_proteome = "host_custom_proteome.fasta"

                #File discordances_after_bmt_simulation = "discordances_after_bmt_simulation.txt"

                File all_discordances = "all_discordances.txt"

                File specific_discordances = tissue+"_specific_discordances.txt"
                File gvhd_specific_discordances = "GvHD_specific_discordances.txt"

                File GvL_dict_nums = "GvL_dict_nums.csv"
                File GvHD_dict_nums = "GvHD_dict_nums.csv"

        }
}


task blastFileGenerationForBlast {
        File specific_discordances

        command {
                # Python file outputs a fa file that is used in the blastOperations task
                # Function that filters the peptide_fasta_id and peptide (pep) from the tissue specific discordances file

                # 2 different specific_discordances are used - mentioned above in the call part 
                python3 /pipeline/blast_preprocessing.py -specific_discordances ${specific_discordances} 
        }

        output {
                File peptides_for_blastp="peptides_for_blast.fa"
        }

        runtime {
                docker : "nidhihookeri/minors-pipeline-2"
                cpu : 16
        }

        meta {
                author: "Nidhi Hookeri"
                email: "nhookeri@broadinstitute.org"
        }
}


task blastOperations {
        File proteome_filename
        String makeblastdb_out
        String blastp_out
        File peptide_file
        String sample

        command {
                # Creates a blast database for the donor sample custom proteome from the task bmt
                # Important to create the database and the run blastp in the same folder as result 
                makeblastdb -in ${proteome_filename} -dbtype prot -input_type fasta -out ${makeblastdb_out}
                # custome proteome db against tissue specific discordances peptide 
                blastp -query ${peptide_file} -db ${makeblastdb_out} -out ${blastp_out} -outfmt "10 sseqid qseqid slen qlen sstart send qstart qend sseq qseq evalue score length pident nident"
        }

        output {
                File phr="${makeblastdb_out}.phr"
                File pin="${makeblastdb_out}.pin"
                File psq="${makeblastdb_out}.psq"
                File blastp_result=sample+"_blastp_peptides_out.csv"
        }

        runtime {
                docker : "biocontainers/blast:v2.2.31_cv2"
                cpu : 16
        }
}


task blastOperationsPostProcessing {
        File specific_discordances
        File input_blastp_result 
        String sample_type
        String sample
        String condition

        command {
                # Assigning column names to the specific_discordances + filter_blast_search 
                # For donor : specific_discordances = tissue_specific_discordances.txt -> donor_discordances_after_blast.txt
                # For host : specific_discordances = donor_discordances_after_blast.txt with genes to remove -> host_discordances_after_blast.txt
                python3 /pipeline/blast_postprocessing.py -specific_discordances ${specific_discordances} -blastp_result ${input_blastp_result} -sample ${sample} -sample_type ${sample_type} -condition ${condition}
        }

        output {
                File discordances_after_blast=sample+"_discordances_after_blast.txt"
                File? genes_with_minors="genes_with_minors.txt"
                File? allgenes_pumas="allgenes_pumas.txt"
                File? GvL_postBlast_dict_nums="GvL_postBlast_dict_nums.csv"
                File? GvHD_postBlast_dict_nums="GvHD_postBlast_dict_nums.csv"
        }

        runtime {
                docker :  "nidhihookeri/minors-pipeline-2"
                cpu : 16
                disks: "local-disk 15 HDD"
        }

        meta {
                author: "Nidhi Hookeri"
                email: "nhookeri@broadinstitute.org"
        }
}



task HLAthena_preprocessing {
        File discordances_after_blast 

        command <<<
                python3 /pipeline/HLAthena_preprocessing.py -discordances_after_blast ${discordances_after_blast}
        >>>

        runtime {
                docker :  "nidhihookeri/minors-pipeline-2"
        }

        output {
                File HLAthena_preprocessed_peptides="HLAthena_peptides.txt"
        }
}



task HLAthena_postprocessing {
        File sample_predictions
        File discordances_after_blast
        String tissue
        String sample_name

        command <<<
                python3 /pipeline/HLAthena_postprocessing.py -hlathena_predictions ${sample_predictions} -discordances ${discordances_after_blast} -tissue_type ${tissue} -sample_name ${sample_name}
        >>>
        
        output {
                File binding_putativeMinorAntigens=sample_name+"_binding_putativeMinorAntigens.txt"
                File? GvL_postHLAthena_dict_nums="GvL_postHLAthena_dict_nums.csv"
                File? GvHD_postHLAthena_dict_nums="GvHD_postHLAthena_dict_nums.txt"
        }

        runtime{
                docker : "nidhihookeri/minors-pipeline-2"
                cpu : 32
        }
}


task merge_summary_csv_files {
    File bcf_dict
    File bmt_gvl_dict
    File bmt_gvhd_dict
    File gvl_postBlast_dict
    File gvhd_postBlast_dict
    File gvl_postHLAthena_dict
    File gvhd_postHLAthena_dict

    String sample_id
    String output_summary_file = sample_id + ".csv"
    String output_summary_file_2 = sample_id + "_format2" + ".csv"

    command <<<
    paste -d, ${bcf_dict} ${bmt_gvl_dict} ${gvl_postBlast_dict} ${gvl_postHLAthena_dict} ${bmt_gvhd_dict} ${gvhd_postBlast_dict} ${gvhd_postHLAthena_dict} >${output_summary_file}

    cat ${bcf_dict} ${bmt_gvl_dict} ${gvl_postBlast_dict} ${gvl_postHLAthena_dict} ${bmt_gvhd_dict} ${gvhd_postBlast_dict} ${gvhd_postHLAthena_dict} > ${output_summary_file_2}
    >>>

    runtime {
      docker : "ubuntu:latest"
    }

    output {
      File output_summary_file_result="${output_summary_file}"
      File output_summary_file_result_2="${output_summary_file_2}"
    }
}