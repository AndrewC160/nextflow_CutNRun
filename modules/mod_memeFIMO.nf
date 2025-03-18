#!/usr/bin/env nextflow

process memeFIMO {
  tag "${sys_idx}"
  cpus 1
  memory '64GB'
  
  publishDir "${params.dir_pool}/${samp_idx}/meme/fimo", mode: 'copy', pattern: "*"
  
  input:
    tuple val(sys_idx), val(samp_idx), val(samp_name), path(fasta_file)
    path motif_file
    val prefix
  
  output:
    path "*"
  
  script:
  """
  fimo --oc ./ ${motif_file} ${fasta_file}
  for file in \$(ls ./*); do
    mv \$file ${prefix}_\$(basename \$file)
  done
  """
}