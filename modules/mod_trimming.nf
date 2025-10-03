#!/usr/bin/env nextflow

process trimming {
  tag "${meta.replicate_name}"
  publishDir "${meta.rep_dir}/qc", mode: 'copy', pattern: "*_trim_report.txt"
  publishDir "${meta.rep_dir}/qc", mode: 'copy', pattern: "*.zip"
  publishDir "${meta.rep_dir}/qc", mode: 'copy', pattern: "*.html"
  
  input:
    tuple val(meta), path(fastq_r1), path(fastq_r2)

  output:
	tuple val(meta), path(fastq_r1_out), path(fastq_r2_out), emit: "trimmed"
	tuple val(meta), path(fastq_r1_out), path(fastq_r2_out), emit: "fastqc"
	tuple path(rpt_fl), path(fastq_r1_fqc), path(fastq_r2_fqc), emit: "report_qc"
	path "*"
  
  script:
  fastq_r1_out = "${meta.replicate_name}_val_1.fq"
  fastq_r2_out = "${meta.replicate_name}_val_2.fq"
  fastq_r1_fqc = "${meta.replicate_name}_val_1_fastqc.zip"
  fastq_r2_fqc = "${meta.replicate_name}_val_2_fastqc.zip"
  rpt_fl = "${meta.replicate_name}_trim_report.txt"
  """
  trim_galore \
    --paired \
    --quality 30 \
    --fastqc \
    --dont_gzip \
    --basename '${meta.replicate_name}' \
    --cores 4 \
    ${fastq_r1} ${fastq_r2}
  
  # Concatentate trimming reports.
  echo -e "==R1==" > ${rpt_fl}
  cat ./*_R1_*.txt >> ${rpt_fl}
  echo -e "==R2==" >> ${rpt_fl}
  cat ./*_R2_*.txt >> ${rpt_fl}
  """
}