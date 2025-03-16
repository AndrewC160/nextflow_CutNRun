#!/usr/bin/env nextflow

process peakCallingBroadPooled{
  tag "${samp_idx}"
  cpus 1
  memory '16GB'
  
  publishDir "${params.dir_pool}/${samp_idx}/qc", mode: 'copy', pattern: "*_report.txt"
  publishDir "${params.dir_pool}/${samp_idx}/qc", mode: 'copy', pattern: "*.tsv"
  publishDir "${params.dir_pool}/${samp_idx}/peaks", mode: 'copy', pattern: "*.broadPeak"
  
  input:
    tuple val(sys_idx), val(samp_idx), val(samp_name), val(proj), val(cell_line), val(epitope), val(condition), path(bams_input), path(bams_ctrl)
    path blacklist_bed
    path seqsize_tsv
  
  output:
    tuple val(sys_idx), val(samp_idx), val(samp_name), val(proj), val(cell_line), val(epitope), val(condition), path("${samp_name}_peaks.broadPeak"), emit: "broadPeaks"
    path "*.txt"
    path "*.tsv"
  
  script:
  rpt_fl = "${samp_name}_broadPeaks_report.txt"
  pks1 = "${samp_name}_all_peaks.broadPeak"
  pks2 = "${samp_name}_peaks.broadPeak"
  rpt_blacklist = "${samp_name}_blacklist_broadPeaks.tsv"
  
  """
  macs3 callpeak \
    --broad \
  	-t ${bams_input} \
  	-c ${bams_ctrl} \
  	-f BAMPE \
  	-g 2.7e9 \
  	-n ${samp_name}'_all'\
  	-q 0.01 \
  	-B \
  	--keep-dup all 2> ${rpt_fl}

  # Remove peaks that overlap blacklisted regions.
  bedtools subtract -A -a ${pks1} -b ${blacklist_bed} > ${pks2}
  
  # Count blacklisted peaks.
  wc -l ${pks1} > ${rpt_blacklist}
  wc -l ${pks2} >> ${rpt_blacklist}
  
  # Remove intermediate files so they don't get published.
  rm ${pks1}
  """
}