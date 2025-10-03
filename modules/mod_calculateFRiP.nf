#!/usr/bin/env nextflow

process calculateFRiP {
  tag "${rep_name}"
  publishDir "${rep_dir}/qc", mode: 'copy'
  
  input:
    tuple val(epitope), val(rep_name), path(bam_file), val(meta), path(bed_file)
    val(suffix)
  
  output:
    path outp_txt, emit: "report_qc"
  
  script:
  rep_name = meta.replicate_name
  rep_dir = meta.rep_dir
  outp_txt = "${rep_name}_${suffix}_FRiP.txt"
  """
  echo -e "total\tin_peaks\tfrip" > ${outp_txt}
  #Count reads.
  reads_total=\$(samtools view -c ${bam_file})
  #reads_total=\$(samtools view --require-flags 3 -c ${bam_file})
  
  #Count reads in peaks.
  reads_peaks=\$(bedtools intersect -a <(cut -f1,2,3 ${bed_file}) -b ${bam_file} -c | awk '{i+=\$4}END{print i}')
  
  #Calculate FRiP.
  frip_val=\$(echo \"scale=5; \$reads_peaks / (\$reads_total)\" | bc)
  echo -e "\$reads_total\t\$reads_peaks\t\$frip_val" >> ${outp_txt}
  """
}