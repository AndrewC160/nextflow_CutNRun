#!/usr/bin/env nextflow

process calculateFRiP {
  tag "${samp_name}"
  publishDir "${params.dir_reps}/${samp_name}/qc", mode: 'copy'
  
  input:
    tuple val(sys_idx), val(samp_idx), val(samp_name), path(bam_file), path(bed_file)
    val(suffix)
  
  output:
    tuple val(sys_idx), val(samp_idx), val(samp_name), path("${samp_name}_${suffix}_FRiP.txt"), emit: "frip"
  
  script:
  outp_txt = "${samp_name}_${suffix}_FRiP.txt"
  """
  echo -e "total\tin_peaks\tfrip" > ${outp_txt}
  
  #Count reads.
  reads_total=\$(samtools view -c ${bam_file})
  
  #Count reads in peaks.
  reads_peaks=\$(bedtools intersect -a <(cut -f1,2,3 ${bed_file}) -b ${bam_file} -c | awk '{i+=\$4}END{print i}')
  
  #Calculate FRiP.
  frip_val=\$(echo \"scale=5; \$reads_peaks / (\$reads_total * 0.5)\" | bc)
  
  echo -e "\$reads_total\t\$reads_peaks\t\$frip_val" >> ${outp_txt}
  """
}