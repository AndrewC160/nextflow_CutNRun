#!/usr/bin/env nextflow

process peakCallingNarrowPooled {
  tag "${samp_name}"
  cpus 1
  memory '16GB'
  
  publishDir "${params.dir_pool}/${samp_name}/qc", mode: 'copy', pattern: "*_report.txt"
  publishDir "${params.dir_pool}/${samp_name}/qc", mode: 'copy', pattern: "*.tsv"
  publishDir "${params.dir_pool}/${samp_name}/peaks", mode: 'copy', pattern: "*.gz"
  publishDir "${params.dir_pool}/${samp_name}/peaks", mode: 'copy', pattern: "*.tbi"
  publishDir "${params.dir_pool}/${samp_name}/peaks", mode: 'copy', pattern: "*.bed"
  publishDir "${params.dir_pool}/${samp_name}/peaks", mode: 'copy', pattern: "*.tsv"
  publishDir "${params.dir_pool}/${samp_name}/peaks", mode: 'copy', pattern: "*.narrowPeak"
  
  input:
    tuple val(proj), val(samp_name), val(cell_line), val(epitope), val(cond), path(bams_input), path(bams_ctrl)
    path blacklist_bed
    path seqsize_tsv
  
  output:
    tuple val(proj), val(samp_name), val(cell_line), val(epitope), val(cond), path("${samp_name}_peaks.narrowPeak"), emit: "narrowPeaks"
    tuple val(proj), val(samp_name), val(cell_line), val(epitope), val(cond), path("${samp_name}_summits.bed"), emit: "summits"
    tuple val(proj), val(samp_name), val(cell_line), val(epitope), val(cond), path("${samp_name}_control_lambda.bdg.gz"), path("${samp_name}_control_lambda.bdg.gz.tbi"), emit: "ctrlBDG"
    tuple val(proj), val(samp_name), val(cell_line), val(epitope), val(cond), path("${samp_name}_treat_pileup.bdg.gz"), path("${samp_name}_treat_pileup.bdg.gz.tbi"), emit: "treatBDG"
    path "*.txt"
    path "*.tsv"
  
  script:
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
  	-t ${bams_input} \
  	-c ${bams_ctrl} \
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