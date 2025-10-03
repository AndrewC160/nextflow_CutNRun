#!/usr/bin/env nextflow

process peakCallingNarrowPooled {
  tag "${pool_name}"
  cpus 1
  memory '16GB'
  
  publishDir "${pool_dir}/qc", mode: 'copy', pattern: "*_report.txt"
  publishDir "${pool_dir}/qc", mode: 'copy', pattern: "*.tsv"
  publishDir "${pool_dir}/peaks", mode: 'copy', pattern: "*.tbi"
  publishDir "${pool_dir}/peaks", mode: 'copy', pattern: "*.gz"
  publishDir "${pool_dir}/peaks", mode: 'copy', pattern: "*.bed"
  publishDir "${pool_dir}/peaks", mode: 'copy', pattern: "*.tsv"
  publishDir "${pool_dir}/peaks", mode: 'copy', pattern: "*.narrowPeak"
  
  input:
	tuple val(meta), path(bams_input), path(bams_ctrl)
    path blacklist_bed
    path seqsize_tsv
  
  output:
	tuple val(meta), path(pks2), emit: "narrowPeaks"
	tuple val(meta), path(sums3), emit: "summits"
	tuple val(meta), path(sums1), emit: "summitsAll"
	tuple val(meta), path(bdg_ctrl3), path(bdg_ctrl3i), emit: "ctrlBDG"
	tuple val(meta), path(bdg_treat3), path(bdg_treat3i), emit: "treatBDG"
    tuple val(meta), path(bams_input), path(bams_ctrl), 
	  path(pks2),
	  path(bdg_treat3), path(bdg_treat3i), 
	  path(bdg_ctrl3), path(bdg_ctrl3i), emit: "ROSE"
	tuple path(pks2), path(sums3), path(bdg_treat3), path(bdg_treat3i), path(bdg_ctrl3), path(bdg_ctrl3i), emit: "report_peaks"
	path rpt_fl, emit: "report_qc"
    path "*.txt"
    path "*.tsv"
  
  script:
  pool_name = meta.pool_name
  pool_dir = meta.pool_dir
  rpt_fl = "${pool_name}_narrowPeaks_report.txt"
  pks1 = "${pool_name}_all_peaks.narrowPeak"
  pks2 = "${pool_name}_peaks.narrowPeak"
  sums1 = "${pool_name}_all_summits.bed"
  sums2 = "${pool_name}_blacklist_summits.bed"
  sums3 = "${pool_name}_summit_regions.bed"
  bdg_ctrl1 = "${pool_name}_all_control_lambda.bdg"
  bdg_ctrl2 = "${pool_name}_control_lambda.bdg"
  bdg_ctrl3 = "${pool_name}_control_lambda.bdg.gz"
  bdg_ctrl3i= "${bdg_ctrl3}.tbi"
  bdg_treat1 = "${pool_name}_all_treat_pileup.bdg"
  bdg_treat2 = "${pool_name}_treat_pileup.bdg"
  bdg_treat3 = "${pool_name}_treat_pileup.bdg.gz"
  bdg_treat3i= "${bdg_treat3}.tbi"
  rpt_blacklist = "${pool_name}_blacklist_narrowPeaks.tsv"
  """
  macs3 callpeak \
  	-t ${bams_input} \
  	-c ${bams_ctrl} \
  	-f BAMPE \
  	-g 2.7e9 \
  	-n ${pool_name}'_all'\
  	-q 0.01 \
  	-B \
  	--call-summits \
  	--keep-dup all 2> ${rpt_fl}

  # Remove peaks that overlap blacklisted regions.
  bedtools subtract -A -a ${pks1} -b ${blacklist_bed} > ${pks2}
  bedtools subtract -A -a ${sums1} -b ${blacklist_bed} > ${sums2}
  
  # Count blacklisted peaks.
  wc -l ${pks1} >  ${rpt_blacklist}
  wc -l ${pks2} >> ${rpt_blacklist}
  
  # Slop summit regions to 500bp windows.
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