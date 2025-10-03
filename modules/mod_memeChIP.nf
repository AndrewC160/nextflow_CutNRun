#!/usr/bin/env nextflow

process memeChIP {
  tag "${meta.pool_name}"
  cpus 1
  memory '64GB'
  
  publishDir "${meta.pool_dir}/meme/memechip/${prefix}", mode: 'copy', pattern: "*"
  
  input:
	tuple val(meta), path(fasta_file)
    path motif_file
    val prefix
  
  output:
    path "*"
  
  script:
  """
  meme-chip --oc ./ -db ${motif_file} ${fasta_file}
  """
}