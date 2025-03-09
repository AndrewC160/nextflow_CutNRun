#!/usr/bin/env nextflow

process memeCENTRIMO {
  tag "${samp_name}"
  cpus 1
  memory '16GB'
  
  publishDir "${params.dir_pool}/${samp_name}/meme/centrimo", mode: 'copy', pattern: "*"
  
  input:
    tuple val(samp_name), val(cell_line), val(epitope), val(condition), path(fasta_file)
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