#!/usr/bin/env nextflow

process memeSEA {
  tag "${meta.pool_name}"
  publishDir "${meta.pool_dir}/meme/sea/${prefix}", mode: 'copy', pattern: "*"
  
  input:
	tuple val(meta), path(fasta_file)
    path motif_file
    val prefix
  
  output:
    path "*"
  
  script:
  """
  sea --p ${fasta_file} --m ${motif_file} --oc ./
  
  # || echo "processed"
  """
}