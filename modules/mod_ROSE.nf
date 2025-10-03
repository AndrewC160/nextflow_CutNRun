#!/usr/bin/env nextflow

process ROSE {
  tag "${meta.pool_name}"
  cpus 6
  memory '16GB'
  publishDir "${meta.pool_dir}/ROSE", mode: 'copy'
  
  input:
	tuple val(meta),
		path(bams_input),
		path(bams_ctrl),
		path(peaks_input),
		path(treat_pileup),
		path(treat_pileup_tbi),
		path(ctrl_pileup),
		path(ctrl_pileup_tbi)
    path(r_markdown)
    path(R_dir)
    path(gene_gtf)
    path(gene_gtf_idx)
  
  output:
    path "*.html"
    path "*.tsv"
    path "*.rds"
    path "*.png"
  
  script:
  rose_file="${meta.pool_name}_ROSE.html"
  """
  Rscript -e 'rmarkdown::render("${r_markdown}",output_format="html_document",
    output_file="${rose_file}",output_dir=".",intermediates_dir=".",quiet=TRUE,
    params=list(
      pool_name="${meta.pool_name}",
      bams_test="${bams_input}",
      bams_ctrl="${bams_ctrl}",
      peak_file="${peaks_input}",
      gene_gtf="${gene_gtf}",
      bdg_treat="${treat_pileup}",
      bdg_ctrl="${ctrl_pileup}",
      R_dir="${R_dir}")
  )'
  """
}