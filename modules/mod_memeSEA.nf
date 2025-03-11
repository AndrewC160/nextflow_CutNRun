#!/usr/bin/env nextflow

process memeSEA {
  tag "${samp_name}"
  publishDir "${params.dir_pool}/${samp_name}/meme/sea", mode: 'copy', pattern: "*"
  
  input:
    tuple val(proj), val(samp_name), val(cell_line), val(epitope), val(condition), path(fasta_file)
    path motif_file
    val prefix
  
  output:
    path "*"
  
  script:
  """
  #NOTE: Error with SEA; appears to run but fails to make HTML. Revisit later.
  sea --p ${fasta_file} --m ${motif_file} --oc ./ || echo "processed"
  
  for file in \$(ls ./*); do
    mv \$file ${prefix}_\$(basename \$file)
  done
  """
}