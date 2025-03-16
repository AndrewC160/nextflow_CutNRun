#!/usr/bin/env nextflow

process getSequences {
  tag "${samp_name}"
  cpus 1
  memory '16GB'
  
  publishDir "${params.dir_pool}/${sys_idx}/peaks", mode: 'copy', pattern: "*.fasta"
  
  input:
    tuple val(sys_idx), val(samp_idx), val(samp_name), val(proj), val(cell_line), val(epitope), val(condition), path(bed_file)
    path genome_fasta
    val prefix
  
  output:
    tuple val(sys_idx), val(samp_idx), val(samp_name), val(proj), val(cell_line), val(epitope), val(condition), path("${samp_name}_${prefix}.fasta"), emit: "seqs"
    path "*"
  
  script:
  out_bam = "${samp_name}_${prefix}.fasta"
  
  """
  bedtools getfasta -name -fi ${genome_fasta} -bed ${bed_file} -fo ${out_bam}
  """
}