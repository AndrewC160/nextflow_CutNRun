#!/usr/bin/env nextflow

process getSequences {
  tag "${meta.pool_name}"
  cpus 1
  memory '16GB'
  
  publishDir "${meta.pool_dir}/peaks", mode: 'copy', pattern: "*.fasta"
  
  input:
	tuple val(meta), path(bed_file)
    path genome_fasta
    val prefix
  
  output:
    tuple val(meta), path(out_fasta), emit: "seqs"
    path "*"
  
  script:
  out_fasta = "${meta.pool_name}_${prefix}.fasta"
  """
  pk_cnt=\$(wc -l <${bed_file})
  if [ \$pk_cnt -eq 1 ]; then
	touch ${out_fasta}
  else
	bedtools getfasta -name -fi ${genome_fasta} -bed ${bed_file} -fo ${out_fasta}
  fi
  """
}