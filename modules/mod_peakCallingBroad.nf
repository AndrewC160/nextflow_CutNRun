#!/usr/bin/env nextflow

process peakCallingBroad {
  tag "${samp_name}"
  cpus 1
  memory '16GB'
  
  publishDir "${out_dir}/${samp_name}/qc", mode: 'copy', pattern: "*_report.txt"
  publishDir "${out_dir}/${samp_name}/qc", mode: 'copy', pattern: "*.tsv"
  publishDir "${out_dir}/${samp_name}/peaks", mode: 'copy', pattern: "*.broadPeak"
  
  input:
    tuple val(samp_name), val(cell_line), val(epitope), val(cond), val(rep), path(bam_file)
    path blacklist_bed
    path seqsize_tsv
    val out_dir
  
  output:
    tuple val(samp_name), val(cell_line), val(epitope), val(cond), val(rep), path("${samp_name}_peaks.broadPeak"), emit: "broadPeaks"
    path "*.txt"
    path "*.tsv"
  
  script:
  bam_input = bam_file
  rpt_fl = "${samp_name}_broadPeaks_report.txt"
  pks1 = "${samp_name}_all_peaks.broadPeak"
  pks2 = "${samp_name}_peaks.broadPeak"
  rpt_blacklist = "${samp_name}_blacklist_broadPeaks.tsv"
  
  """
  macs3 callpeak \
    --broad \
  	-t ${bam_input} \
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