#!/usr/bin/env nextflow

process getSummitRegions {
  tag "${samp_idx}"
  publishDir "${out_dir}/${samp_name}/peaks", mode: 'copy', pattern: "*.bed"
  
  input:
    tuple val(samp_idx), val(samp_name), path(summit_file)
  
  output:
    tuple val(samp_idx), path("${samp_name}_summit_regions.tsv"), emit: "spikeTable"
  
  script:
  tsv_spike = "${samp_name}_spike.tsv"
  """
  Rscript ${params.dir_R}/combine_spike_tsvs.R ${spike_files} > "${samp_name}_spike.tsv"
  """
}