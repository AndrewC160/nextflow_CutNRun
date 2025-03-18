#!/usr/bin/env nextflow

process peakCallingNarrow {
  tag "${samp_name}"
  cpus 1
  memory '16GB'
  
  publishDir "${out_dir}/${samp_name}/qc", mode: 'copy', pattern: "*_report.txt"
  publishDir "${out_dir}/${samp_name}/qc", mode: 'copy', pattern: "*.tsv"
  publishDir "${out_dir}/${samp_name}/peaks", mode: 'copy', pattern: "*.gz"
  publishDir "${out_dir}/${samp_name}/peaks", mode: 'copy', pattern: "*.tbi"
  publishDir "${out_dir}/${samp_name}/peaks", mode: 'copy', pattern: "*.bed"
  publishDir "${out_dir}/${samp_name}/peaks", mode: 'copy', pattern: "*.tsv"
  publishDir "${out_dir}/${samp_name}/peaks", mode: 'copy', pattern: "*.narrowPeak"
  
  input:
    tuple val(sys_idx), val(samp_idx), val(cond_name), val(samp_name), val(epitope), path(bam_file)
    path blacklist_bed
    path seqsize_tsv
    val out_dir
  
  output:
    tuple val(sys_idx), val(samp_idx), val(cond_name), val(samp_name), val(epitope), path("${samp_name}_peaks.narrowPeak"), emit: "narrowPeaks"
    tuple val(sys_idx), val(samp_idx), val(cond_name), val(samp_name), val(epitope), path("${samp_name}_summits.bed"), emit: "summits"
    tuple val(sys_idx), val(samp_idx), val(cond_name), val(samp_name), val(epitope), path("${samp_name}_control_lambda.bdg.gz"), path("${samp_name}_control_lambda.bdg.gz.tbi"), emit: "ctrlBDG"
    tuple val(sys_idx), val(samp_idx), val(cond_name), val(samp_name), val(epitope), path("${samp_name}_treat_pileup.bdg.gz"), path("${samp_name}_treat_pileup.bdg.gz.tbi"), emit: "treatBDG"
    path "*.txt"
    path "*.tsv"
  
  script:
  bam_input = bam_file
  rpt_fl = "${samp_name}_narrowPeaks_report.txt"
  pks1 = "${samp_name}_all_peaks.narrowPeak"
  pks2 = "${samp_name}_peaks.narrowPeak"
  sums1 = "${samp_name}_all_summits.bed"
  sums2 = "${samp_name}_blacklist_summits.bed"
  sums3 = "${samp_name}_summits.bed"
  bdg_ctrl1 = "${samp_name}_all_control_lambda.bdg"
  bdg_ctrl2 = "${samp_name}_control_lambda.bdg"
  bdg_ctrl3 = "${samp_name}_control_lambda.bdg.gz"
  bdg_treat1 = "${samp_name}_all_treat_pileup.bdg"
  bdg_treat2 = "${samp_name}_treat_pileup.bdg"
  bdg_treat3 = "${samp_name}_treat_pileup.bdg.gz"
  rpt_blacklist = "${samp_name}_blacklist_narrowPeaks.tsv"
  """
  macs3 callpeak \
  	-t ${bam_input} \
  	-f BAMPE \
  	-g 2.7e9 \
  	-n ${samp_name}'_all'\
  	-q 0.01 \
  	-B \
  	--call-summits \
  	--keep-dup all 2> ${rpt_fl}

  # Remove peaks that overlap blacklisted regions.
  bedtools subtract -A -a ${pks1} -b ${blacklist_bed} > ${pks2}
  bedtools subtract -A -a ${sums1} -b ${blacklist_bed} > ${sums2}
  
  # Count blacklisted peaks.
  wc -l ${pks1} > ${rpt_blacklist}
  wc -l ${pks2} >> ${rpt_blacklist}
  
  # Slop summit regions to 201bp windows.
  bedtools slop -i ${sums2} -g ${seqsize_tsv} -b 100 > ${sums3} || true
  
  # Rename, BGZip, and index bedgraph files.
  mv ${bdg_ctrl1} ${bdg_ctrl2}
  mv ${bdg_treat1} ${bdg_treat2}
  
  bgzip ${bdg_ctrl2} &&
  bgzip ${bdg_treat2} &&
  wait
  
  tabix -p 'bed' ${bdg_ctrl3} &&
  tabix -p 'bed' ${bdg_treat3} &&
  wait
  
  # Remove intermediate files so they don't get published.
  rm ${pks1} ${sums1} ${sums2}
  """
}