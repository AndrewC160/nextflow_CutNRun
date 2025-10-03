#!/usr/bin/env nextflow

process memeCENTRIMO {
  tag "${meta.pool_name}"
  cpus 1
  memory '16GB'
  
  publishDir "${meta.pool_dir}/meme/centrimo/${prefix}", mode: 'copy', pattern: "*"
  
  input:
	tuple val(meta), path(fasta_file)
    path motif_file
    val prefix
  
  output:
    path "*"
  
  script:
  """
  centrimo --oc ./ ${fasta_file} ${motif_file}
  #|| echo "processed"
  """
}