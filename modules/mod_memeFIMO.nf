#!/usr/bin/env nextflow

process memeFIMO {
  tag "${samp_name}"
  cpus 1
  memory '64GB'
  
  publishDir "${params.dir_pool}/${samp_name}/meme/fimo", mode: 'copy', pattern: "*"
  
  input:
    tuple val(samp_name), val(cell_line), val(epitope), val(condition), path(fasta_file)
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



//74/dd2b4e