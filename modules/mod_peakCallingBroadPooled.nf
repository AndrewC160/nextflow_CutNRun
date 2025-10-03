#!/usr/bin/env nextflow

process peakCallingBroadPooled{
  tag "${pool_name}"
  cpus 1
  memory '16GB'
  
  publishDir "${pool_dir}/qc", mode: 'copy', pattern: "*_report.txt"
  publishDir "${pool_dir}/qc", mode: 'copy', pattern: "*.tsv"
  publishDir "${pool_dir}/peaks", mode: 'copy', pattern: "*.broadPeak"
  
  input:
	tuple val(meta), path(bams_input), path(bams_ctrl)
    path blacklist_bed
    path seqsize_tsv
  
  output:
	tuple val(meta), path(pks2), emit: "broadPeaks"
	path pks2, emit: "report_peaks"
	path rpt_fl, emit: "report_qc"
    path "*.txt"
    path "*.tsv"
  
  script:
  pool_name = meta.pool_name
  pool_dir = meta.pool_dir
  rpt_fl = "${pool_name}_broadPeaks_report.txt"
  pks1 = "${pool_name}_all_peaks.broadPeak"
  pks2 = "${pool_name}_peaks.broadPeak"
  rpt_blacklist = "${pool_name}_blacklist_broadPeaks.tsv"
  """
  macs3 callpeak \
    --broad \
  	-t ${bams_input} \
  	-c ${bams_ctrl} \
  	-f BAMPE \
  	-g 2.7e9 \
  	-n ${pool_name}'_all'\
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