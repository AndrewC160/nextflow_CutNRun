#!/usr/bin/env nextflow

process combineSpikes {
  tag "${cond}"
  publishDir "${params.dir_pool}/${cond}", mode: 'copy', pattern: "*.tsv"
  
  input:
    tuple val(cond), path(spike_files)
  
  output:
    tuple val(cond), path("${cond}_spike.tsv"), emit: "spikeTable"
  
  script:
  tsv_spike = "${cond}_spike.tsv"
  """
  Rscript ${params.dir_R}/combine_spike_tsvs.R ${spike_files} > "${cond}_spike.tsv"
  """
}