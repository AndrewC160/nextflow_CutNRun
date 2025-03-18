#!/usr/bin/env nextflow

process poolReport {
  tag "${pool_name}"
  publishDir "${params.dir_pool}/${pool_name}", mode: 'copy', pattern: "*.html"
  
  input:
    path(r_markdown)
    tuple val(pool_name), 
          path(fqc_trim1),
          path(fqc_trim2),
          path(fqc_filt),
          path(tsv_spikes),
          path(npks_rep),
          path(bdgs_rep),
          path(bdgs_rep_tbi),
          path(bdgs_rep_ctrl),
          path(bdgs_rep_ctrl_tbi),
          path(bpks_rep),
          path(frip),
          path(npks_pool),
          path(bpks_pool),
          path(bdgs_pool),
          path(bdgs_pool_tbi),
          path(bdgs_pool_ctrl),
          path(bdgs_pool_ctrl_tbi),
          path(sums_pool)
    val(ctrl_epitope)
    path(sample_table)
    val(output_dir)
    path(R_dir)
    path(gene_gtf)
    path(gene_gtf_idx)
  
  output:
    path "${pool_name}_report.html"
  
  script:
  report_file="${pool_name}_report.html"
  """
  Rscript -e 'rmarkdown::render("${r_markdown}",output_format="html_document",
    output_file="${report_file}",output_dir=".",intermediates_dir=".",quiet=TRUE,
    params=list(
      pool_name="${pool_name}",
      ctrl_epitope="${ctrl_epitope}",
      fastqc_trim1="${fqc_trim1}",
      fastqc_trim2="${fqc_trim2}",
      fastqc_filt="${fqc_filt}",
      tsv_spike="${tsv_spikes}",
      npks_reps="${npks_rep}",
      bdgs_reps="${bdgs_rep}",
      bdgs_reps_ctrl="${bdgs_rep_ctrl}",
      bpks_reps="${bpks_rep}",
      tsvs_frip="${frip}",
      npks_pool="${npks_pool}",
      sums_pool="${sums_pool}",
      bdgs_pool="${bdgs_pool}",
      bdgs_pool_ctrl="${bdgs_pool_ctrl}",
      bpks_pool="${bpks_pool}",
      sample_table="${sample_table}",
      output_dir="${output_dir}",
      gene_gtf="${gene_gtf}",
      R_dir="${R_dir}"
    )
  )'
  """
}