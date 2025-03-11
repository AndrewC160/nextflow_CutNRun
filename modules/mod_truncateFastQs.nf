#!/usr/bin/env nextflow

process truncateFastQs{
  tag "${samp_name}"
  cpus 6
  
  input:
  tuple val(proj), val(samp_name), val(cell_line), val(epitope), val(cond), val(rep), path(fastq_r1), path(fastq_r2)
  val fastq_reads
  
  output:
  tuple val(proj), val(samp_name), val(cell_line), val(epitope), val(cond), val(rep), path(fastq_r1_out), path(fastq_r2_out)
  
  script:
  fastq_r1_out = "${fastq_r1.getSimpleName()}_trunc.fastq.gz"
  fastq_r2_out = "${fastq_r2.getSimpleName()}_trunc.fastq.gz"
  """
  zcat ${fastq_r1} | head -n \$((${fastq_reads}*4)) | gzip -c > ${fastq_r1_out} &
  zcat ${fastq_r2} | head -n \$((${fastq_reads}*4)) | gzip -c > ${fastq_r2_out} &
  wait
  """
}