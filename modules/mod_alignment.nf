#!/usr/bin/env nextflow

process alignment {
  tag "${samp_name}"
  cpus 8
  memory '32.GB'
  
  publishDir "${params.dir_reps}/${samp_name}/qc", mode: 'copy', pattern: "*.txt"
  
  input:
    tuple val(sys_idx), val(samp_idx), val(cond_name), val(proj), val(samp_name), val(cell_line), val(epitope), val(cond), val(rep), path(fastq_r1), path(fastq_r2)
    path bowtie2_idx
  
  output:
    tuple val(sys_idx), val(samp_idx), val(cond_name), val(proj), val(samp_name), val(cell_line), val(epitope), val(cond), val(rep), path("*.bam"), emit: "aligned"
    path "*"
  
  script:
  gen_nm = bowtie2_idx.baseName
  out_bam = "${samp_name}_${gen_nm}_unfiltered.bam"
  rpt_fl = "${samp_name}_${gen_nm}_alignment_report.txt"
  
  """
  bowtie2 \
    -x ${bowtie2_idx}/${gen_nm} \
    -1 ${fastq_r1} \
    -2 ${fastq_r2} \
    --dovetail \
    --reorder \
    --threads 8 \
    2> ${rpt_fl} | \
    samtools view -bS - > ${out_bam}
  """
}