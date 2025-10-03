#!/usr/bin/env nextflow

process poolReport {
  tag "${proj_name}"
  publishDir "${params.dir_out}/${proj_name}", mode: 'copy', pattern: "*.html"
  publishDir "${params.dir_out}/${proj_name}", mode: 'copy', pattern: "*.tsv"
  
  input:
	val proj_name
	path reps_tsv
	path pool_tsv
	path sample_table
	val ctrl_epitope
	val output_dir
	path R_dir
	path seqsizes, stageAs: "seqsizes.tsv"
	path gene_gtf, stageAs: "genes.gtf"
	path gene_gtf_idx, stageAs: "genes.gtf.tbi"
	path r_markdown
	path trim_rep, stageAs: "reps/qc/*"
	path alg_prime, stageAs: "reps/align/*"
	path alg_spike, stageAs: "reps/align/*"
	path bam_best, stageAs: "reps/align/*"
	path flt_prime, stageAs: "reps/qc/*"
	path flt_prime_alg, stageAs: "reps/align/*"
	path npk_reps, stageAs: "reps/qc/*"
	path npk_reps_pks, stageAs: "reps/peaks/*"
	path bpk_reps, stageAs: "reps/qc/*"
	path bpk_reps_pks, stageAs: "reps/peaks/*"
	path npk_pool, stageAs: "pool/qc/*"
	path npk_pool_pks, stageAs: "pool/peaks/*"
	path bpk_pool, stageAs: "pool/qc/*"
	path bpk_pool_pks, stageAs: "pool/peaks/*"
	path frip_reps, stageAs: "reps/qc/*"
	path spikes_pool, stageAs: "pool/align/*"
	//path fimo, stageAs: "pool/fimo/*"
	//path fimo_npks, stageAs: "pool/fimo_npks/*"
	//path sea, stageAs: "pool/sea/*"
	//path sea_npks, stageAs: "pool/sea_npks/*"
	//path centrimo, stageAs: "pool/centrimo/*"
	//path rose, stageAs: "pool/rose/*"

  output:
    path report_file
	path "*.tsv"

  script:
  report_file="${proj_name}_report.html"
  
  """
  Rscript -e 'rmarkdown::render("${r_markdown}",output_format="html_document",
	output_file="${report_file}",output_dir=".",intermediates_dir=".",quiet=TRUE,
	params=list(
		ctrl_epitope="${ctrl_epitope}",
		out_dir="${output_dir}"
	)
  )'
  """
}