version 1.0

# =============================================================================
# MASPHL viral NGS end-to-end workflow (SC2 / MPXV / WNV)
# -----------------------------------------------------------------------------
# Phase 2: load from BaseSpace + analyze with TheiaCoV. (Phase 3 adds QC
# PASS/FAIL and the result/LIMS files.)
#
# Launch in Terra with "Run workflow with inputs defined by file paths".
#
#   1. basespace_load       check inputs, check run, register samples, fetch reads
#                           (workflows/wf_basespace_load.wdl, tested in Phase 1)
#   2. TheiaCoV v4.0.0      per sample, imported unmodified from Theiagen's PHB
#   3. MPXV only            if the all-clades Nextclade call is I, Ia or Ib, rerun
#                           TheiaCoV with the clade I reference, primer bed and
#                           Nextclade dataset; the rerun's results are kept
#   4. write_results        every TheiaCoV output -> a data table column, same
#                           names as Terra's default output mapping
#
# Extra TheiaCoV settings (anything not listed in the inputs below) can still be
# set in Terra's Inputs tab under the task name "theiacov_illumina_pe" (and
# "theiacov_illumina_pe_clade_i" for the MPXV rerun), as in today's configs.
# =============================================================================

import "wf_basespace_load.wdl" as loader
import "../tasks/task_terra_table.wdl" as terra_table
# TheiaCoV_Illumina_PE_PHB v4.0.0 (Theiagen, Public Health Bioinformatics), unmodified.
import "https://raw.githubusercontent.com/theiagen/public_health_bioinformatics/v4.0.0/workflows/theiacov/wf_theiacov_illumina_pe.wdl" as phb_theiacov

workflow masphl_viral_ngs {
  meta {
    description: "MASPHL SC2/MPXV/WNV: load a BaseSpace run, run TheiaCoV (with MPXV clade I rerun), write results to the data table."
  }
  parameter_meta {
    basespace_run_id: "Number from the run's BaseSpace URL: basespace.illumina.com/run/<number>/details"
    btb_batch_id: "SC2 only: 5-digit BtB batch ID."
    only_these_samples: "Optional. Comma-separated sample IDs to load instead of the whole project. Controls are not added automatically."
    extra_negative_control_ids: "Optional. Comma-separated IDs of controls whose names don't start with NEG or NTC."
    organism: "TheiaCoV organism: sars-cov-2 (SC2), MPXV or WNV."
    mpxv_clade_i_reference_genome: "MPXV workspace: reference used to rerun clade I samples."
  }

  input {
    # ---- per run ---------------------------------------------------------
    String basespace_run_id
    String? btb_batch_id
    String? only_these_samples
    String? extra_negative_control_ids

    # ---- workspace defaults: pathogen, BaseSpace, Terra ---------------------
    String pathogen
    String? basespace_project
    String access_token
    String api_server = "https://api.basespace.illumina.com"
    String terra_project
    String terra_workspace

    # ---- workspace defaults: TheiaCoV (same values as today's TheiaCoV config)
    String organism
    File? reference_genome
    File? primer_bed
    File? reference_gff
    File? reference_gene_locations_bed
    Int? genome_length
    String? nextclade_dataset_name
    String? nextclade_dataset_tag

    # ---- MPXV workspace only: clade I rerun --------------------------------
    File? mpxv_clade_i_reference_genome
    File? mpxv_clade_i_primer_bed
    File? mpxv_clade_i_reference_gff
    File? mpxv_clade_i_reference_gene_locations_bed
    Int? mpxv_clade_i_genome_length
    String mpxv_clade_i_nextclade_dataset_name = "nextstrain/mpox/clade-i"
    String? mpxv_clade_i_nextclade_dataset_tag

    # ---- rarely changed ----------------------------------------------------
    String table_name = "sample"
    Boolean write_to_table = true
    String? sample_set_name
    String control_id_regex = "^(NEG|NTC)"
    Boolean allow_incomplete_run = false
    Boolean clear_old_columns_on_rerun = true
    String masphl_workflow_version = "phase2-dev"
  }

  # ---- configuration problems, reported by the loader's quick input check ---
  Boolean organism_ok =
    if pathogen == "SC2" then (organism == "sars-cov-2" || organism == "SARS-CoV-2")
    else if pathogen == "MPXV" then (organism == "MPXV" || organism == "mpox" || organism == "Mpox" ||
                                     organism == "monkeypox" || organism == "Monkeypox virus")
    else if pathogen == "WNV" then (organism == "WNV" || organism == "wnv" || organism == "West Nile virus")
    else true   # an unknown pathogen is reported by the loader itself
  Boolean is_mpxv = pathogen == "MPXV"
  Array[String] config_problems = flatten([
    if organism_ok then [] else
      ["organism '" + organism + "' does not match pathogen '" + pathogen + "'. Check the workspace defaults."],
    if is_mpxv && !defined(nextclade_dataset_name) then
      ["nextclade_dataset_name must be set for MPXV (e.g. nextstrain/mpox/all-clades); without it TheiaCoV uses a dataset that never reports clade I."] else [],
    if is_mpxv && !defined(mpxv_clade_i_reference_genome) then
      ["mpxv_clade_i_reference_genome must be set for MPXV runs (used to rerun clade I samples)."] else [],
    if is_mpxv && !defined(mpxv_clade_i_primer_bed) then
      ["mpxv_clade_i_primer_bed must be set for MPXV runs (used to rerun clade I samples)."] else []
  ])

  # ===========================================================================
  # 1. Load from BaseSpace
  # ===========================================================================
  call loader.basespace_load as load {
    input:
      basespace_run_id = basespace_run_id,
      btb_batch_id = btb_batch_id,
      only_these_samples = only_these_samples,
      extra_negative_control_ids = extra_negative_control_ids,
      pathogen = pathogen,
      basespace_project = basespace_project,
      access_token = access_token,
      api_server = api_server,
      terra_project = terra_project,
      terra_workspace = terra_workspace,
      table_name = table_name,
      write_to_table = write_to_table,
      sample_set_name = sample_set_name,
      control_id_regex = control_id_regex,
      allow_incomplete_run = allow_incomplete_run,
      clear_old_columns_on_rerun = clear_old_columns_on_rerun,
      loader_version = masphl_workflow_version,
      other_input_problems = config_problems
  }

  # ===========================================================================
  # 2-3. TheiaCoV per sample, with the MPXV clade I rerun
  # ===========================================================================
  scatter (i in range(length(load.loaded_sample_ids))) {
    String sample_id = load.loaded_sample_ids[i]

    call phb_theiacov.theiacov_illumina_pe as theiacov_illumina_pe {
      input:
        samplename = sample_id,
        read1 = load.read1[i],
        read2 = load.read2[i],
        organism = organism,
        reference_genome = reference_genome,
        primer_bed = primer_bed,
        reference_gff = reference_gff,
        reference_gene_locations_bed = reference_gene_locations_bed,
        genome_length = genome_length,
        nextclade_dataset_name = nextclade_dataset_name,
        nextclade_dataset_tag = nextclade_dataset_tag
    }

    # TheiaCoV writes Nextclade's "clade" column as-is ("NA" when unassigned).
    # Clade I = I, Ia or Ib (also accepts a "/..." suffix, e.g. "Ib/sh2023").
    String first_clade = select_first([theiacov_illumina_pe.nextclade_clade, "NA"])
    Boolean rerun_clade_i = is_mpxv && sub(first_clade, "^I[ab]?(/.*)?$", "CLADE_I") == "CLADE_I"

    if (rerun_clade_i) {
      call phb_theiacov.theiacov_illumina_pe as theiacov_illumina_pe_clade_i {
        input:
          samplename = sample_id,
          read1 = load.read1[i],
          read2 = load.read2[i],
          organism = organism,
          reference_genome = mpxv_clade_i_reference_genome,
          primer_bed = mpxv_clade_i_primer_bed,
          reference_gff = mpxv_clade_i_reference_gff,
          reference_gene_locations_bed = mpxv_clade_i_reference_gene_locations_bed,
          genome_length = mpxv_clade_i_genome_length,
          nextclade_dataset_name = mpxv_clade_i_nextclade_dataset_name,
          nextclade_dataset_tag = mpxv_clade_i_nextclade_dataset_tag
      }
    }

    String analysis_route = if !is_mpxv then "standard"
                            else if rerun_clade_i then "mpxv_clade_i_rerun"
                            else "mpxv_all_clades"

    # BEGIN GENERATED theiacov values (tools/gen_theiacov_outputs.py; do not edit by hand)
    Array[String] theiacov_values = [
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.theiacov_illumina_pe_version}" else "~{theiacov_illumina_pe.theiacov_illumina_pe_version}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.theiacov_illumina_pe_analysis_date}" else "~{theiacov_illumina_pe.theiacov_illumina_pe_analysis_date}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.seq_platform}" else "~{theiacov_illumina_pe.seq_platform}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.read_screen_raw}" else "~{theiacov_illumina_pe.read_screen_raw}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.read_screen_raw_tsv}" else "~{theiacov_illumina_pe.read_screen_raw_tsv}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.read_screen_clean}" else "~{theiacov_illumina_pe.read_screen_clean}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.read_screen_clean_tsv}" else "~{theiacov_illumina_pe.read_screen_clean_tsv}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_num_reads_raw1}" else "~{theiacov_illumina_pe.fastq_scan_num_reads_raw1}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_num_reads_raw2}" else "~{theiacov_illumina_pe.fastq_scan_num_reads_raw2}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_num_reads_raw_pairs}" else "~{theiacov_illumina_pe.fastq_scan_num_reads_raw_pairs}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_r1_mean_readlength_raw}" else "~{theiacov_illumina_pe.fastq_scan_r1_mean_readlength_raw}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_r2_mean_readlength_raw}" else "~{theiacov_illumina_pe.fastq_scan_r2_mean_readlength_raw}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_r1_mean_q_raw}" else "~{theiacov_illumina_pe.fastq_scan_r1_mean_q_raw}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_r2_mean_q_raw}" else "~{theiacov_illumina_pe.fastq_scan_r2_mean_q_raw}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_version}" else "~{theiacov_illumina_pe.fastq_scan_version}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_docker}" else "~{theiacov_illumina_pe.fastq_scan_docker}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_raw1_json}" else "~{theiacov_illumina_pe.fastq_scan_raw1_json}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_raw2_json}" else "~{theiacov_illumina_pe.fastq_scan_raw2_json}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_num_reads_clean1}" else "~{theiacov_illumina_pe.fastq_scan_num_reads_clean1}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_num_reads_clean2}" else "~{theiacov_illumina_pe.fastq_scan_num_reads_clean2}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_num_reads_clean_pairs}" else "~{theiacov_illumina_pe.fastq_scan_num_reads_clean_pairs}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_r1_mean_readlength_clean}" else "~{theiacov_illumina_pe.fastq_scan_r1_mean_readlength_clean}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_r2_mean_readlength_clean}" else "~{theiacov_illumina_pe.fastq_scan_r2_mean_readlength_clean}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_r1_mean_q_clean}" else "~{theiacov_illumina_pe.fastq_scan_r1_mean_q_clean}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_r2_mean_q_clean}" else "~{theiacov_illumina_pe.fastq_scan_r2_mean_q_clean}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_clean1_json}" else "~{theiacov_illumina_pe.fastq_scan_clean1_json}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastq_scan_clean2_json}" else "~{theiacov_illumina_pe.fastq_scan_clean2_json}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastqc_num_reads_raw1}" else "~{theiacov_illumina_pe.fastqc_num_reads_raw1}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastqc_num_reads_raw2}" else "~{theiacov_illumina_pe.fastqc_num_reads_raw2}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastqc_num_reads_raw_pairs}" else "~{theiacov_illumina_pe.fastqc_num_reads_raw_pairs}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastqc_num_reads_clean1}" else "~{theiacov_illumina_pe.fastqc_num_reads_clean1}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastqc_num_reads_clean2}" else "~{theiacov_illumina_pe.fastqc_num_reads_clean2}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastqc_num_reads_clean_pairs}" else "~{theiacov_illumina_pe.fastqc_num_reads_clean_pairs}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastqc_raw1_html}" else "~{theiacov_illumina_pe.fastqc_raw1_html}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastqc_raw2_html}" else "~{theiacov_illumina_pe.fastqc_raw2_html}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastqc_clean1_html}" else "~{theiacov_illumina_pe.fastqc_clean1_html}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastqc_clean2_html}" else "~{theiacov_illumina_pe.fastqc_clean2_html}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastqc_version}" else "~{theiacov_illumina_pe.fastqc_version}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastqc_docker}" else "~{theiacov_illumina_pe.fastqc_docker}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.trimmomatic_version}" else "~{theiacov_illumina_pe.trimmomatic_version}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.trimmomatic_docker}" else "~{theiacov_illumina_pe.trimmomatic_docker}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastp_version}" else "~{theiacov_illumina_pe.fastp_version}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.fastp_html_report}" else "~{theiacov_illumina_pe.fastp_html_report}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.read1_clean}" else "~{theiacov_illumina_pe.read1_clean}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.read2_clean}" else "~{theiacov_illumina_pe.read2_clean}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.bbduk_docker}" else "~{theiacov_illumina_pe.bbduk_docker}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.read1_dehosted}" else "~{theiacov_illumina_pe.read1_dehosted}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.read2_dehosted}" else "~{theiacov_illumina_pe.read2_dehosted}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.kraken_version}" else "~{theiacov_illumina_pe.kraken_version}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.kraken_human}" else "~{theiacov_illumina_pe.kraken_human}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.kraken_sc2}" else "~{theiacov_illumina_pe.kraken_sc2}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.kraken_target_organism}" else "~{theiacov_illumina_pe.kraken_target_organism}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.kraken_target_organism_name}" else "~{theiacov_illumina_pe.kraken_target_organism_name}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.kraken_report}" else "~{theiacov_illumina_pe.kraken_report}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.kraken_human_dehosted}" else "~{theiacov_illumina_pe.kraken_human_dehosted}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.kraken_sc2_dehosted}" else "~{theiacov_illumina_pe.kraken_sc2_dehosted}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.kraken_target_organism_dehosted}" else "~{theiacov_illumina_pe.kraken_target_organism_dehosted}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.kraken_report_dehosted}" else "~{theiacov_illumina_pe.kraken_report_dehosted}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.bwa_version}" else "~{theiacov_illumina_pe.bwa_version}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.samtools_version}" else "~{theiacov_illumina_pe.samtools_version}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.read1_aligned}" else "~{theiacov_illumina_pe.read1_aligned}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.read2_aligned}" else "~{theiacov_illumina_pe.read2_aligned}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.aligned_bam}" else "~{theiacov_illumina_pe.aligned_bam}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.aligned_bai}" else "~{theiacov_illumina_pe.aligned_bai}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.read1_unaligned}" else "~{theiacov_illumina_pe.read1_unaligned}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.read2_unaligned}" else "~{theiacov_illumina_pe.read2_unaligned}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.sorted_bam_unaligned}" else "~{theiacov_illumina_pe.sorted_bam_unaligned}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.sorted_bam_unaligned_bai}" else "~{theiacov_illumina_pe.sorted_bam_unaligned_bai}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.primer_trimmed_read_percent}" else "~{theiacov_illumina_pe.primer_trimmed_read_percent}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.ivar_version_primtrim}" else "~{theiacov_illumina_pe.ivar_version_primtrim}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.samtools_version_primtrim}" else "~{theiacov_illumina_pe.samtools_version_primtrim}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.primer_bed_name}" else "~{theiacov_illumina_pe.primer_bed_name}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.ivar_tsv}" else "~{theiacov_illumina_pe.ivar_tsv}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.ivar_vcf}" else "~{theiacov_illumina_pe.ivar_vcf}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.ivar_variant_proportion_intermediate}" else "~{theiacov_illumina_pe.ivar_variant_proportion_intermediate}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.ivar_variant_version}" else "~{theiacov_illumina_pe.ivar_variant_version}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.samtools_version_consensus}" else "~{theiacov_illumina_pe.samtools_version_consensus}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.assembly_method}" else "~{theiacov_illumina_pe.assembly_method}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.assembly_fasta}" else "~{theiacov_illumina_pe.assembly_fasta}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.ivar_version_consensus}" else "~{theiacov_illumina_pe.ivar_version_consensus}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.consensus_n_variant_min_depth}" else "~{theiacov_illumina_pe.consensus_n_variant_min_depth}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.consensus_stats}" else "~{theiacov_illumina_pe.consensus_stats}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.consensus_flagstat}" else "~{theiacov_illumina_pe.consensus_flagstat}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.meanbaseq_trim}" else "~{theiacov_illumina_pe.meanbaseq_trim}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.meanmapq_trim}" else "~{theiacov_illumina_pe.meanmapq_trim}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.assembly_mean_coverage}" else "~{theiacov_illumina_pe.assembly_mean_coverage}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.samtools_version_stats}" else "~{theiacov_illumina_pe.samtools_version_stats}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.number_N}" else "~{theiacov_illumina_pe.number_N}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.assembly_length_unambiguous}" else "~{theiacov_illumina_pe.assembly_length_unambiguous}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.number_Degenerate}" else "~{theiacov_illumina_pe.number_Degenerate}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.number_Total}" else "~{theiacov_illumina_pe.number_Total}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.percent_reference_coverage}" else "~{theiacov_illumina_pe.percent_reference_coverage}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.sc2_s_gene_mean_coverage}" else "~{theiacov_illumina_pe.sc2_s_gene_mean_coverage}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.sc2_s_gene_percent_coverage}" else "~{theiacov_illumina_pe.sc2_s_gene_percent_coverage}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.est_percent_gene_coverage_tsv}" else "~{theiacov_illumina_pe.est_percent_gene_coverage_tsv}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.pango_lineage}" else "~{theiacov_illumina_pe.pango_lineage}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.pango_lineage_expanded}" else "~{theiacov_illumina_pe.pango_lineage_expanded}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.pangolin_conflicts}" else "~{theiacov_illumina_pe.pangolin_conflicts}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.pangolin_notes}" else "~{theiacov_illumina_pe.pangolin_notes}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.pangolin_assignment_version}" else "~{theiacov_illumina_pe.pangolin_assignment_version}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.pango_lineage_report}" else "~{theiacov_illumina_pe.pango_lineage_report}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.pangolin_docker}" else "~{theiacov_illumina_pe.pangolin_docker}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.pangolin_versions}" else "~{theiacov_illumina_pe.pangolin_versions}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_version}" else "~{theiacov_illumina_pe.nextclade_version}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_docker}" else "~{theiacov_illumina_pe.nextclade_docker}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_json}" else "~{theiacov_illumina_pe.nextclade_json}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.auspice_json}" else "~{theiacov_illumina_pe.auspice_json}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_tsv}" else "~{theiacov_illumina_pe.nextclade_tsv}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_ds_tag}" else "~{theiacov_illumina_pe.nextclade_ds_tag}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_aa_subs}" else "~{theiacov_illumina_pe.nextclade_aa_subs}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_aa_dels}" else "~{theiacov_illumina_pe.nextclade_aa_dels}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_clade}" else "~{theiacov_illumina_pe.nextclade_clade}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_lineage}" else "~{theiacov_illumina_pe.nextclade_lineage}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_qc}" else "~{theiacov_illumina_pe.nextclade_qc}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_json_flu_h5n1}" else "~{theiacov_illumina_pe.nextclade_json_flu_h5n1}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.auspice_json_flu_h5n1}" else "~{theiacov_illumina_pe.auspice_json_flu_h5n1}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_tsv_flu_h5n1}" else "~{theiacov_illumina_pe.nextclade_tsv_flu_h5n1}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_aa_subs_flu_h5n1}" else "~{theiacov_illumina_pe.nextclade_aa_subs_flu_h5n1}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_aa_dels_flu_h5n1}" else "~{theiacov_illumina_pe.nextclade_aa_dels_flu_h5n1}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_clade_flu_h5n1}" else "~{theiacov_illumina_pe.nextclade_clade_flu_h5n1}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_qc_flu_h5n1}" else "~{theiacov_illumina_pe.nextclade_qc_flu_h5n1}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_json_flu_ha}" else "~{theiacov_illumina_pe.nextclade_json_flu_ha}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.auspice_json_flu_ha}" else "~{theiacov_illumina_pe.auspice_json_flu_ha}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_tsv_flu_ha}" else "~{theiacov_illumina_pe.nextclade_tsv_flu_ha}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_ds_tag_flu_ha}" else "~{theiacov_illumina_pe.nextclade_ds_tag_flu_ha}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_aa_subs_flu_ha}" else "~{theiacov_illumina_pe.nextclade_aa_subs_flu_ha}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_aa_dels_flu_ha}" else "~{theiacov_illumina_pe.nextclade_aa_dels_flu_ha}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_clade_flu_ha}" else "~{theiacov_illumina_pe.nextclade_clade_flu_ha}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_qc_flu_ha}" else "~{theiacov_illumina_pe.nextclade_qc_flu_ha}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_json_flu_na}" else "~{theiacov_illumina_pe.nextclade_json_flu_na}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.auspice_json_flu_na}" else "~{theiacov_illumina_pe.auspice_json_flu_na}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_tsv_flu_na}" else "~{theiacov_illumina_pe.nextclade_tsv_flu_na}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_ds_tag_flu_na}" else "~{theiacov_illumina_pe.nextclade_ds_tag_flu_na}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_aa_subs_flu_na}" else "~{theiacov_illumina_pe.nextclade_aa_subs_flu_na}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_aa_dels_flu_na}" else "~{theiacov_illumina_pe.nextclade_aa_dels_flu_na}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_clade_flu_na}" else "~{theiacov_illumina_pe.nextclade_clade_flu_na}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.nextclade_qc_flu_na}" else "~{theiacov_illumina_pe.nextclade_qc_flu_na}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.vadr_alerts_list}" else "~{theiacov_illumina_pe.vadr_alerts_list}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.vadr_feature_tbl_pass}" else "~{theiacov_illumina_pe.vadr_feature_tbl_pass}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.vadr_feature_tbl_fail}" else "~{theiacov_illumina_pe.vadr_feature_tbl_fail}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.vadr_classification_summary_file}" else "~{theiacov_illumina_pe.vadr_classification_summary_file}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.vadr_all_outputs_tar_gz}" else "~{theiacov_illumina_pe.vadr_all_outputs_tar_gz}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.vadr_num_alerts}" else "~{theiacov_illumina_pe.vadr_num_alerts}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.vadr_docker}" else "~{theiacov_illumina_pe.vadr_docker}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.vadr_fastas_zip_archive}" else "~{theiacov_illumina_pe.vadr_fastas_zip_archive}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.irma_version}" else "~{theiacov_illumina_pe.irma_version}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.irma_docker}" else "~{theiacov_illumina_pe.irma_docker}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.irma_type}" else "~{theiacov_illumina_pe.irma_type}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.irma_subtype}" else "~{theiacov_illumina_pe.irma_subtype}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.irma_subtype_notes}" else "~{theiacov_illumina_pe.irma_subtype_notes}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.irma_assembly_fasta_concatenated}" else "~{theiacov_illumina_pe.irma_assembly_fasta_concatenated}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.irma_ha_segment_fasta}" else "~{theiacov_illumina_pe.irma_ha_segment_fasta}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.irma_na_segment_fasta}" else "~{theiacov_illumina_pe.irma_na_segment_fasta}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.irma_pa_segment_fasta}" else "~{theiacov_illumina_pe.irma_pa_segment_fasta}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.irma_pb1_segment_fasta}" else "~{theiacov_illumina_pe.irma_pb1_segment_fasta}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.irma_pb2_segment_fasta}" else "~{theiacov_illumina_pe.irma_pb2_segment_fasta}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.irma_mp_segment_fasta}" else "~{theiacov_illumina_pe.irma_mp_segment_fasta}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.irma_np_segment_fasta}" else "~{theiacov_illumina_pe.irma_np_segment_fasta}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.irma_ns_segment_fasta}" else "~{theiacov_illumina_pe.irma_ns_segment_fasta}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.irma_qc_summary_tsv}" else "~{theiacov_illumina_pe.irma_qc_summary_tsv}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.irma_all_snvs_tsv}" else "~{theiacov_illumina_pe.irma_all_snvs_tsv}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.irma_all_insertions_tsv}" else "~{theiacov_illumina_pe.irma_all_insertions_tsv}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.irma_all_deletions_tsv}" else "~{theiacov_illumina_pe.irma_all_deletions_tsv}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.genoflu_version}" else "~{theiacov_illumina_pe.genoflu_version}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.genoflu_genotype}" else "~{theiacov_illumina_pe.genoflu_genotype}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.genoflu_all_segments}" else "~{theiacov_illumina_pe.genoflu_all_segments}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.genoflu_output_tsv}" else "~{theiacov_illumina_pe.genoflu_output_tsv}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.abricate_flu_type}" else "~{theiacov_illumina_pe.abricate_flu_type}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.abricate_flu_subtype}" else "~{theiacov_illumina_pe.abricate_flu_subtype}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.abricate_flu_results}" else "~{theiacov_illumina_pe.abricate_flu_results}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.abricate_flu_database}" else "~{theiacov_illumina_pe.abricate_flu_database}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.abricate_flu_version}" else "~{theiacov_illumina_pe.abricate_flu_version}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.flu_A_315675_resistance}" else "~{theiacov_illumina_pe.flu_A_315675_resistance}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.flu_amantadine_resistance}" else "~{theiacov_illumina_pe.flu_amantadine_resistance}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.flu_compound_367_resistance}" else "~{theiacov_illumina_pe.flu_compound_367_resistance}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.flu_favipiravir_resistance}" else "~{theiacov_illumina_pe.flu_favipiravir_resistance}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.flu_fludase_resistance}" else "~{theiacov_illumina_pe.flu_fludase_resistance}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.flu_L_742_001_resistance}" else "~{theiacov_illumina_pe.flu_L_742_001_resistance}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.flu_laninamivir_resistance}" else "~{theiacov_illumina_pe.flu_laninamivir_resistance}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.flu_peramivir_resistance}" else "~{theiacov_illumina_pe.flu_peramivir_resistance}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.flu_pimodivir_resistance}" else "~{theiacov_illumina_pe.flu_pimodivir_resistance}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.flu_rimantadine_resistance}" else "~{theiacov_illumina_pe.flu_rimantadine_resistance}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.flu_oseltamivir_resistance}" else "~{theiacov_illumina_pe.flu_oseltamivir_resistance}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.flu_xofluza_resistance}" else "~{theiacov_illumina_pe.flu_xofluza_resistance}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.flu_zanamivir_resistance}" else "~{theiacov_illumina_pe.flu_zanamivir_resistance}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.quasitools_version}" else "~{theiacov_illumina_pe.quasitools_version}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.quasitools_date}" else "~{theiacov_illumina_pe.quasitools_date}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.quasitools_coverage_file}" else "~{theiacov_illumina_pe.quasitools_coverage_file}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.quasitools_dr_report}" else "~{theiacov_illumina_pe.quasitools_dr_report}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.quasitools_hydra_vcf}" else "~{theiacov_illumina_pe.quasitools_hydra_vcf}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.quasitools_mutations_report}" else "~{theiacov_illumina_pe.quasitools_mutations_report}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.qc_check}" else "~{theiacov_illumina_pe.qc_check}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.qc_standard}" else "~{theiacov_illumina_pe.qc_standard}",
      if (rerun_clade_i) then "~{theiacov_illumina_pe_clade_i.percentage_mapped_reads}" else "~{theiacov_illumina_pe.percentage_mapped_reads}"
    ]
    # END GENERATED theiacov values

    Array[String] result_row = flatten([
      [sample_id, analysis_route,
       first_clade, "~{theiacov_illumina_pe.nextclade_lineage}",
       masphl_workflow_version],
      theiacov_values
    ])
  }

  # BEGIN GENERATED theiacov columns (tools/gen_theiacov_outputs.py; do not edit by hand)
  Array[String] theiacov_columns = [
    "theiacov_illumina_pe_version",
    "theiacov_illumina_pe_analysis_date",
    "seq_platform",
    "read_screen_raw",
    "read_screen_raw_tsv",
    "read_screen_clean",
    "read_screen_clean_tsv",
    "fastq_scan_num_reads_raw1",
    "fastq_scan_num_reads_raw2",
    "fastq_scan_num_reads_raw_pairs",
    "fastq_scan_r1_mean_readlength_raw",
    "fastq_scan_r2_mean_readlength_raw",
    "fastq_scan_r1_mean_q_raw",
    "fastq_scan_r2_mean_q_raw",
    "fastq_scan_version",
    "fastq_scan_docker",
    "fastq_scan_raw1_json",
    "fastq_scan_raw2_json",
    "fastq_scan_num_reads_clean1",
    "fastq_scan_num_reads_clean2",
    "fastq_scan_num_reads_clean_pairs",
    "fastq_scan_r1_mean_readlength_clean",
    "fastq_scan_r2_mean_readlength_clean",
    "fastq_scan_r1_mean_q_clean",
    "fastq_scan_r2_mean_q_clean",
    "fastq_scan_clean1_json",
    "fastq_scan_clean2_json",
    "fastqc_num_reads_raw1",
    "fastqc_num_reads_raw2",
    "fastqc_num_reads_raw_pairs",
    "fastqc_num_reads_clean1",
    "fastqc_num_reads_clean2",
    "fastqc_num_reads_clean_pairs",
    "fastqc_raw1_html",
    "fastqc_raw2_html",
    "fastqc_clean1_html",
    "fastqc_clean2_html",
    "fastqc_version",
    "fastqc_docker",
    "trimmomatic_version",
    "trimmomatic_docker",
    "fastp_version",
    "fastp_html_report",
    "read1_clean",
    "read2_clean",
    "bbduk_docker",
    "read1_dehosted",
    "read2_dehosted",
    "kraken_version",
    "kraken_human",
    "kraken_sc2",
    "kraken_target_organism",
    "kraken_target_organism_name",
    "kraken_report",
    "kraken_human_dehosted",
    "kraken_sc2_dehosted",
    "kraken_target_organism_dehosted",
    "kraken_report_dehosted",
    "bwa_version",
    "samtools_version",
    "read1_aligned",
    "read2_aligned",
    "aligned_bam",
    "aligned_bai",
    "read1_unaligned",
    "read2_unaligned",
    "sorted_bam_unaligned",
    "sorted_bam_unaligned_bai",
    "primer_trimmed_read_percent",
    "ivar_version_primtrim",
    "samtools_version_primtrim",
    "primer_bed_name",
    "ivar_tsv",
    "ivar_vcf",
    "ivar_variant_proportion_intermediate",
    "ivar_variant_version",
    "samtools_version_consensus",
    "assembly_method",
    "assembly_fasta",
    "ivar_version_consensus",
    "consensus_n_variant_min_depth",
    "consensus_stats",
    "consensus_flagstat",
    "meanbaseq_trim",
    "meanmapq_trim",
    "assembly_mean_coverage",
    "samtools_version_stats",
    "number_N",
    "assembly_length_unambiguous",
    "number_Degenerate",
    "number_Total",
    "percent_reference_coverage",
    "sc2_s_gene_mean_coverage",
    "sc2_s_gene_percent_coverage",
    "est_percent_gene_coverage_tsv",
    "pango_lineage",
    "pango_lineage_expanded",
    "pangolin_conflicts",
    "pangolin_notes",
    "pangolin_assignment_version",
    "pango_lineage_report",
    "pangolin_docker",
    "pangolin_versions",
    "nextclade_version",
    "nextclade_docker",
    "nextclade_json",
    "auspice_json",
    "nextclade_tsv",
    "nextclade_ds_tag",
    "nextclade_aa_subs",
    "nextclade_aa_dels",
    "nextclade_clade",
    "nextclade_lineage",
    "nextclade_qc",
    "nextclade_json_flu_h5n1",
    "auspice_json_flu_h5n1",
    "nextclade_tsv_flu_h5n1",
    "nextclade_aa_subs_flu_h5n1",
    "nextclade_aa_dels_flu_h5n1",
    "nextclade_clade_flu_h5n1",
    "nextclade_qc_flu_h5n1",
    "nextclade_json_flu_ha",
    "auspice_json_flu_ha",
    "nextclade_tsv_flu_ha",
    "nextclade_ds_tag_flu_ha",
    "nextclade_aa_subs_flu_ha",
    "nextclade_aa_dels_flu_ha",
    "nextclade_clade_flu_ha",
    "nextclade_qc_flu_ha",
    "nextclade_json_flu_na",
    "auspice_json_flu_na",
    "nextclade_tsv_flu_na",
    "nextclade_ds_tag_flu_na",
    "nextclade_aa_subs_flu_na",
    "nextclade_aa_dels_flu_na",
    "nextclade_clade_flu_na",
    "nextclade_qc_flu_na",
    "vadr_alerts_list",
    "vadr_feature_tbl_pass",
    "vadr_feature_tbl_fail",
    "vadr_classification_summary_file",
    "vadr_all_outputs_tar_gz",
    "vadr_num_alerts",
    "vadr_docker",
    "vadr_fastas_zip_archive",
    "irma_version",
    "irma_docker",
    "irma_type",
    "irma_subtype",
    "irma_subtype_notes",
    "irma_assembly_fasta_concatenated",
    "irma_ha_segment_fasta",
    "irma_na_segment_fasta",
    "irma_pa_segment_fasta",
    "irma_pb1_segment_fasta",
    "irma_pb2_segment_fasta",
    "irma_mp_segment_fasta",
    "irma_np_segment_fasta",
    "irma_ns_segment_fasta",
    "irma_qc_summary_tsv",
    "irma_all_snvs_tsv",
    "irma_all_insertions_tsv",
    "irma_all_deletions_tsv",
    "genoflu_version",
    "genoflu_genotype",
    "genoflu_all_segments",
    "genoflu_output_tsv",
    "abricate_flu_type",
    "abricate_flu_subtype",
    "abricate_flu_results",
    "abricate_flu_database",
    "abricate_flu_version",
    "flu_A_315675_resistance",
    "flu_amantadine_resistance",
    "flu_compound_367_resistance",
    "flu_favipiravir_resistance",
    "flu_fludase_resistance",
    "flu_L_742_001_resistance",
    "flu_laninamivir_resistance",
    "flu_peramivir_resistance",
    "flu_pimodivir_resistance",
    "flu_rimantadine_resistance",
    "flu_oseltamivir_resistance",
    "flu_xofluza_resistance",
    "flu_zanamivir_resistance",
    "quasitools_version",
    "quasitools_date",
    "quasitools_coverage_file",
    "quasitools_dr_report",
    "quasitools_hydra_vcf",
    "quasitools_mutations_report",
    "qc_check",
    "qc_standard",
    "percentage_mapped_reads"
  ]
  # END GENERATED theiacov columns

  # ===========================================================================
  # 4. Write TheiaCoV results to the data table
  # ===========================================================================
  call terra_table.make_result_rows {
    input:
      table_name = table_name,
      columns = flatten([
        ["analysis_route", "nextclade_clade_all_clades", "nextclade_lineage_all_clades",
         "masphl_workflow_version"],
        theiacov_columns
      ]),
      rows = result_row
  }

  if (write_to_table) {
    call terra_table.update_terra_table as write_results {
      input:
        terra_project = terra_project,
        terra_workspace = terra_workspace,
        table_name = table_name,
        rows_tsv = make_result_rows.rows_tsv,
        check_previous_run = false,
        wait_for = load.read_paths_update_log    # after the loader's writes
    }
  }

  output {
    # loader
    File check_summary = load.check_summary
    File check_report = load.check_report
    String sample_set = load.sample_set
    File? replaced_rows = load.replaced_rows
    # analysis
    Array[String] analyzed_sample_ids = sample_id
    Array[String] analysis_routes = analysis_route
    File theiacov_results = make_result_rows.rows_tsv
    File? results_update_log = write_results.update_log
  }
}
