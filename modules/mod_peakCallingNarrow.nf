#!/usr/bin/env nextflow

process peakCallingNarrow {
  tag "${rep_name}"
  cpus 1
  memory '16GB'
  
  publishDir "${rep_dir}/qc", mode: 'copy', pattern: "*_report.txt"
  publishDir "${rep_dir}/qc", mode: 'copy', pattern: "*.tsv"
  publishDir "${rep_dir}/peaks", mode: 'copy', pattern: "*.gz"
  publishDir "${rep_dir}/peaks", mode: 'copy', pattern: "*.tbi"
  publishDir "${rep_dir}/peaks", mode: 'copy', pattern: "*.bed"
  publishDir "${rep_dir}/peaks", mode: 'copy', pattern: "*.tsv"
  publishDir "${rep_dir}/peaks", mode: 'copy', pattern: "*.narrowPeak"
  
  input:
    tuple val(meta), val(bam_file)
    path blacklist_bed
    path seqsize_tsv
  
  output:
    tuple val(meta), path(pks2), emit: "narrowPeaks"
    tuple val(meta), path(sums3), emit: "summits"
    tuple val(meta), path(sums1), emit: "summitsAll"
    tuple val(meta), path(bdg_ctrl3), path(bdg_ctrl3i), emit: "ctrlBDG"
    tuple val(meta), path(bdg_treat3), path(bdg_treat3i), emit: "treatBDG"
	tuple path(pks2), path(sums3), path(bdg_ctrl3), path(bdg_ctrl3i), path(bdg_treat3), path(bdg_treat3i), emit: "report_peaks"
	path rpt_fl, emit: "report_qc"
    path "*.txt"
    path "*.tsv"
  
  script:
  rep_name= meta.replicate_name
  rep_dir = meta.rep_dir
  rpt_fl = "${rep_name}_narrowPeaks_report.txt"
  pks1 = "${rep_name}_all_peaks.narrowPeak"
  pks2 = "${rep_name}_peaks.narrowPeak"
  sums1 = "${rep_name}_all_summits.bed"
  sums2 = "${rep_name}_blacklist_summits.bed"
  sums3 = "${rep_name}_summit_regions.bed"
  bdg_ctrl1 = "${rep_name}_all_control_lambda.bdg"
  bdg_ctrl2 = "${rep_name}_control_lambda.bdg"
  bdg_ctrl3 = "${rep_name}_control_lambda.bdg.gz"
  bdg_ctrl3i ="${bdg_ctrl3}.tbi"
  bdg_treat1 = "${rep_name}_all_treat_pileup.bdg"
  bdg_treat2 = "${rep_name}_treat_pileup.bdg"
  bdg_treat3 = "${rep_name}_treat_pileup.bdg.gz"
  bdg_treat3i ="${bdg_treat3}.tbi"
  rpt_blacklist = "${rep_name}_blacklist_narrowPeaks.tsv"
  """
  macs3 callpeak \
  	-t ${bam_file} \
  	-f BAMPE \
  	-g 2.7e9 \
  	-n ${rep_name}'_all'\
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
  
  # Slop summit regions to 500bp windows.
  # bedtools slop -i ${sums2} -g ${seqsize_tsv} -b 250 > ${sums3} || true
  Rscript ${params.dir_R}/merge_summits.R ${sums2} ${sums3} 500
  
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
  rm ${pks1} ${sums2}
  """
}