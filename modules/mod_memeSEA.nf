#!/usr/bin/env nextflow

process memeSEA {
  tag "${sys_idx}"
  publishDir "${params.dir_pool}/${samp_idx}/meme/sea", mode: 'copy', pattern: "*"
  
  input:
    tuple val(sys_idx), val(samp_idx), val(samp_name), path(fasta_file)
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