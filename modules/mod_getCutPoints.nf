#!/usr/bin/env nextflow

process getCutPoints {
  tag "${samp_name}"
  publishDir "${out_dir}/${cond_name}/peaks/cuts", mode: 'copy', pattern: "*"
  
  input:
    tuple val(cond_name), path(bam_file), path(bed_file)
    val out_dir
    val prefix
  
  output:
    path "${samp_name}_${prefix}_cuts.bed", emit: "cuts"
  
  script:
  samp_name = "${bam_file.getSimpleName()}"
  tmp_bed = "${bam_file}.tmp"
  out_bed = "${samp_name}_${prefix}_cuts.bed"
  """
  bedtools intersect -wa -a ${bam_file} -b ${bed_file} | \
  	bedtools bamtobed -i stdin > ${tmp_bed}
  
  cat \
    <(awk '{OFS=\"\t\"}{print \$1,\$2,\$2+1,\$4,\$5,\$6}' ${tmp_bed}) \
    <(awk '{OFS=\"\t\"}{print \$1,\$3-1,\$3,\$4,\$5,\$6}' ${tmp_bed}) > ${out_bed}
  """
}