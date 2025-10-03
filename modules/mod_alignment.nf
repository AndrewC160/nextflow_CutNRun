#!/usr/bin/env nextflow

process alignment {
  tag "${meta.replicate_name}"
  cpus 8
  memory '32.GB'
  
  publishDir "${meta.rep_dir}/qc", mode: 'copy', pattern: "*.txt"
  
  input:
    tuple val(meta), path(fastq_r1), path(fastq_r2)
    path bowtie2_idx
  
  output:
    tuple val(meta), path(out_bam), emit: "aligned"
	path rpt_fl, emit: "report_qc"
    path "*"
  
  script:
  gen_nm = bowtie2_idx.baseName
  out_bam = "${meta.replicate_name}_${gen_nm}_unfiltered.bam"
  rpt_fl = "${meta.replicate_name}_${gen_nm}_alignment_report.txt"
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