#!/usr/bin/env nextflow

process readFiltering {
  tag "${replicate_name}"
  cpus 8
  memory '32.GB'
  
  publishDir "${meta.rep_dir}/qc", mode: 'copy', pattern: "*.txt"
  publishDir "${meta.rep_dir}/qc", mode: 'copy', pattern: "*.zip"
  publishDir "${meta.rep_dir}/qc", mode: 'copy', pattern: "*.html"
  publishDir "${meta.rep_dir}/align", mode: 'copy', pattern: "*.bam"
  publishDir "${meta.rep_dir}/align", mode: 'copy', pattern: "*.bai"
  
  input:
    tuple val(meta), val(genome), path(bam_file)
    path blacklist_bed
    
  output:
    tuple val(meta), path(bam_shrt), emit: "filtered"
	tuple val(meta), path(bam_long), emit: "filter_long"
	tuple val(meta), path(fastqc_zip), emit: "fastqc"
	path fastqc_zip, emit: "report_qc"
	path bam_shrt, emit: "report_align"
	path "*.zip"
    path "*.txt"
    path "*.html"
    
  script:
  replicate_name = meta.replicate_name
  pfx = "${replicate_name}_${genome}"
  bam_long = "${pfx}_long.bam"
  bam_shrt = "${pfx}_short.bam"
  rpt_dups = "${pfx}_dup_filter.txt"
  rpt_blist= "${pfx}_blacklist_reads.txt"
  fastqc_zip = "${pfx}_short_fastqc.zip"  
  """
  # Sort by coordinate, put temp files in work directory to facilitate larger files.
  samtools sort -@ 8 -O 'bam' -T 'srt_tmp' ${bam_file} > sorted.bam
  
  # Add @RG line to header so that MarkDuplicates doesn't explode.
  samtools addreplacerg -w -r '@RG\tID:RG1\tSM:${replicate_name}' -o withRG.bam sorted.bam
  
  # Remove duplicates.
  picard MarkDuplicates \
    -I withRG.bam \
    -O deduped.bam \
    -M ${rpt_dups} \
    --ASSUME_SORTED true \
    --REMOVE_DUPLICATES true \
    --VERBOSITY ERROR \
    --QUIET true
    
  # Remove reads in blacklist regions.
  # Note: We don't really care about Sac3 blacklist sites, but when intersected
  # with hg38 blacklist sites this creates an error. Just pipe these into a sink
  # and move on.
  bedtools intersect -v -a withRG.bam -b ${blacklist_bed} > blacklisted.bam 2>/dev/null
  samtools view -c deduped.bam > ${rpt_blist}
  samtools view -c blacklisted.bam >> ${rpt_blist}
  
  # Filter long and short reads into separate files.
  samtools view -h blacklisted.bam | awk 'length(\$10) <= 120 || \$1 ~ /^@/' | samtools view -bS - > ${bam_shrt} & \
    samtools view -h blacklisted.bam | awk 'length(\$10) > 120 || \$1 ~ /^@/' | samtools view -bS - > ${bam_long} & \
    wait

  # Index.
  samtools index ${bam_long} & 
    samtools index ${bam_shrt} &
    wait

  # FastQC
  fastqc ${bam_long} ${bam_shrt}
  """
}