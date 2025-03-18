#!/usr/bin/env nextflow

process memeCENTRIMO {
  tag "${sys_idx}"
  cpus 1
  memory '16GB'
  
  publishDir "${params.dir_pool}/${samp_idx}/meme/centrimo", mode: 'copy', pattern: "*"
  
  input:
    tuple val(sys_idx), val(samp_idx), val(samp_name), path(fasta_file)
    path motif_file
    val prefix
  
  output:
    path "*"
  
  script:
  """
  centrimo --oc ./ ${fasta_file} ${motif_file} || echo "processed"
  for file in \$(ls ./*); do
    mv \$file ${prefix}_\$(basename \$file)
  done
  """
}