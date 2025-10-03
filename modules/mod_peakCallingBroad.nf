#!/usr/bin/env nextflow

process peakCallingBroad {
  tag "${rep_name}"
  cpus 1
  memory '16GB'
  
  publishDir "${rep_dir}/qc", mode: 'copy', pattern: "*_report.txt"
  publishDir "${rep_dir}/qc", mode: 'copy', pattern: "*.tsv"
  publishDir "${rep_dir}/peaks", mode: 'copy', pattern: "*.broadPeak"
  
  input:
    tuple val(meta), val(bam_file)
    path blacklist_bed
    path seqsize_tsv
  
  output:
    tuple val(meta), path(pks2), emit: "broadPeaks"
	path pks2, emit: "report_peaks"
	path rpt_fl, emit: "report_qc"
    path "*.txt"
    path "*.tsv"
  
  script:
  rep_name = meta.replicate_name
  rep_dir = meta.rep_dir
  rpt_fl = "${rep_name}_broadPeaks_report.txt"
  pks1 = "${rep_name}_all_peaks.broadPeak"
  pks2 = "${rep_name}_peaks.broadPeak"
  rpt_blacklist = "${rep_name}_blacklist_broadPeaks.tsv"
  """
  macs3 callpeak \
    --broad \
  	-t ${bam_file} \
  	-f BAMPE \
  	-g 2.7e9 \
  	-n ${rep_name}'_all'\
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