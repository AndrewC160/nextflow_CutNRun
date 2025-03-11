#!/usr/bin/env nextflow

process readFiltering {
  tag "${samp_name}"
  cpus 8
  memory '32.GB'
  
  publishDir "${params.dir_reps}/${samp_name}/qc", mode: 'copy', pattern: "*.txt"
  publishDir "${params.dir_reps}/${samp_name}/qc", mode: 'copy', pattern: "*.zip"
  publishDir "${params.dir_reps}/${samp_name}/qc", mode: 'copy', pattern: "*.html"
  publishDir "${params.dir_reps}/${samp_name}/align", mode: 'copy', pattern: "*.bam"
  publishDir "${params.dir_reps}/${samp_name}/align", mode: 'copy', pattern: "*.bai"
  
  input:
    tuple val(proj), val(samp_name), val(cell_line), val(epitope), val(cond), val(rep), path(bam_file)
    val(gen_nm)
    path(blacklist_bed)
    
  output:
    tuple val(proj), val(samp_name), val(cell_line), val(epitope), val(cond), val(rep), path("${samp_name}_${gen_nm}.bam"), emit: "filtered"
    tuple val(proj), val(samp_name), val(cell_line), val(epitope), val(cond), val(rep), path("${samp_name}_${gen_nm}_long.bam"), emit: "filtered_long"
    path "*.txt"
    path "*.zip"
    path "*.html"
    
  script:
  bam1 = "${bam_file}.sorted"
  bam2 = "${bam_file}.withRG"
  bam3 = "${bam_file}.deduped"
  bam4 = "${bam_file}.blacklisted"
  bam_long = "${samp_name}_${gen_nm}_long.bam"
  bam_shrt = "${samp_name}_${gen_nm}_short.bam"
  rpt_dups = "${samp_name}_${gen_nm}_dup_filter.txt"
  rpt_blist= "${samp_name}_${gen_nm}_blacklist_reads.txt"
  
  """
  # Sort by coordinate, put temp files in work directory to facilitate larger files.
  samtools sort -@ 8 -O 'bam' -T 'srt_tmp' ${bam_file} > ${bam1}
  
  # Add @RG line to header so that MarkDuplicates doesn't explode.
  samtools addreplacerg -w -r '@RG\tID:RG1\tSM:${samp_name}\tPL:Illumina\tLB:Library.fa' -o ${bam2} ${bam1}
  
  # Remove duplicates.
  picard MarkDuplicates \
    -I ${bam2} \
    -O ${bam3} \
    -M ${rpt_dups} \
    -ASSUME_SORTED true \
    -REMOVE_DUPLICATES true \
    --QUIET true
    
  # Remove reads in blacklist regions.
  # Note: We don't really care about Sac3 blacklist sites, but when intersected
  # with hg38 blacklist sites this creates an error. Just pipe these into a sink
  # and move on.
  bedtools intersect -v -a ${bam3} -b ${blacklist_bed} > ${bam4} 2>/dev/null
  samtools view -c ${bam3} > ${rpt_blist}
  samtools view -c ${bam4} >> ${rpt_blist}
  
  # Filter long and short reads into separate files.
  samtools view -h ${bam4} | awk 'length(\$10) <= 120 || \$1 ~ /^@/' | samtools view -bS - > ${bam_shrt} & \
    samtools view -h ${bam4} | awk 'length(\$10) > 120 || \$1 ~ /^@/' | samtools view -bS - > ${bam_long} & \
    wait
  
  # Filter long and short reads into separate files.
  samtools view -h ${bam4} | awk 'length(\$10) <= 120 || \$1 ~ /^@/' | samtools view -bS - > ${bam_shrt} &
    samtools view -h ${bam4} | awk 'length(\$10) > 120 || \$1 ~ /^@/' | samtools view -bS - > ${bam_long} &
    wait

  # Index.
  samtools index ${bam_long} & 
    samtools index ${bam_shrt} &
    wait

  # FastQC
  fastqc ${bam_long} ${bam_shrt}
  """
}