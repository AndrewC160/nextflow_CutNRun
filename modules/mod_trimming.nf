#!/usr/bin/env nextflow

process trimming {
  tag "${samp_name}"
  publishDir "${params.dir_reps}/${samp_name}/qc", mode: 'copy', pattern: "*_trim_report.txt"
  publishDir "${params.dir_reps}/${samp_name}/qc", mode: 'copy', pattern: "*.zip"
  publishDir "${params.dir_reps}/${samp_name}/qc", mode: 'copy', pattern: "*.html"
  
  input:
    tuple val(sys_idx), val(samp_idx), val(cond_name), val(proj), val(samp_name), val(cell_line), val(epitope), val(cond), val(rep), path(fastq_r1), path(fastq_r2)
  
  output:
    tuple val(sys_idx), val(samp_idx), val(cond_name), val(proj), val(samp_name), val(cell_line), val(epitope), val(cond), val(rep), path("${samp_name}_val_1.fq"), path("${samp_name}_val_2.fq"), emit: "trimmed"
    tuple val(sys_idx), val(samp_idx), val(cond_name), val(proj), val(samp_name), val(cell_line), val(epitope), val(cond), val(rep), path("${samp_name}_val_1_fastqc.zip"), path("${samp_name}_val_2_fastqc.zip"), emit: "fastqc"
    path "*"
  
  script:
  rpt_fl = "${samp_name}_trim_report.txt"
  """
  trim_galore \
    --paired \
    --quality 30 \
    --fastqc \
    --dont_gzip \
    --basename '${samp_name}' \
    --cores 4 \
    ${fastq_r1} ${fastq_r2}
  
  # Concatentate trimming reports.
  echo -e "==R1==" > ${rpt_fl}
  cat ./*_R1_*.txt >> ${rpt_fl}
  echo -e "==R2==" >> ${rpt_fl}
  cat ./*_R2_*.txt >> ${rpt_fl}
  """
}