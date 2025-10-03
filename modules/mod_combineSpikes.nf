#!/usr/bin/env nextflow

process combineSpikes {
  tag "${pool_name}"
  publishDir "${pool_dir}/align", mode: 'copy', pattern: "*.tsv"
  
  input:
    tuple val(epitope), val(pool_name), val(pool_dir), path(spike_files)
  
  output:
    tuple val(pool_name), path(tsv_spike), emit: "spikeTable"
	path tsv_spike, emit: "report_align"
  
  script:
  tsv_spike = "${pool_name}_spike.tsv"
  """
  Rscript ${params.dir_R}/combine_spike_tsvs.R ${spike_files} > "${tsv_spike}"
  """
}