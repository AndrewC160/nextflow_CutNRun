#!/usr/bin/env nextflow

process combineSpikes {
  tag "${samp_idx}"
  publishDir "${params.dir_pool}/${samp_idx}", mode: 'copy', pattern: "*.tsv"
  
  input:
    tuple val(samp_idx), val(samp_name), path(spike_files)
  
  output:
    tuple val(samp_idx), path("${samp_name}_spike.tsv"), emit: "spikeTable"
  
  script:
  tsv_spike = "${samp_name}_spike.tsv"
  """
  Rscript ${params.dir_R}/combine_spike_tsvs.R ${spike_files} > "${samp_name}_spike.tsv"
  """
}