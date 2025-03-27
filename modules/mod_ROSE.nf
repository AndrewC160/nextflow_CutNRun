#!/usr/bin/env nextflow

process ROSE {
  tag "${samp_idx}"
  cpus 6
  memory '16GB'
  publishDir "${params.dir_pool}/${samp_idx}/ROSE", mode: 'copy'
  
  input:
    tuple val(sys_idx), val(samp_idx), val(samp_name),
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
  rose_file="${samp_idx}_ROSE.html"
  """
  Rscript -e 'rmarkdown::render("${r_markdown}",output_format="html_document",
    output_file="${rose_file}",output_dir=".",intermediates_dir=".",quiet=TRUE,
    params=list(
      pool_name="${samp_idx}",
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