#!/usr/bin/env nextflow

process bamBest {
  tag "${samp_name}"
  publishDir "${params.dir_reps}/${samp_name}/align", mode: 'copy', pattern: "*.tsv"
  
  input:
    tuple val(samp_name), val(cell_line), val(epitope), val(cond), val(rep), path(bam_hg38_in), path(bam_sac3_in)
  
  output:
    tuple val(samp_name), val(cell_line), val(epitope), val(cond), val(rep), path("${samp_name}_hg38.bam"), emit: "hg38"
    tuple val(samp_name), val(cell_line), val(epitope), val(cond), val(rep), path("${samp_name}_sac3.bam"), emit: "sac3"
    tuple val(samp_name), val(cell_line), val(epitope), val(cond), val(rep), path("${samp_name}_spike.tsv"), emit: "spike"
  
  script:
  bam_hg38_out = "${samp_name}_hg38.bam"
  bam_sac3_out = "${samp_name}_sac3.bam"
  tsv_spike = "${samp_name}_spike.tsv"
  """
  ngsutilsj bam-best \
    --stats ${tsv_spike} \
    --unsorted \
    ${bam_hg38_in} \
    ${bam_sac3_in} \
    -- \
    ${bam_hg38_out} \
    ${bam_sac3_out}
  """
}