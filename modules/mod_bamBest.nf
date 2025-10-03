#!/usr/bin/env nextflow

process bamBest {
  tag "${rep_id}"
  publishDir "${meta.rep_dir}/align", mode: 'copy', pattern: "*.tsv"
  
  input:
	tuple val(rep_id), val(meta), path(bam_primary), val(meta_spike), path(bam_spike)
  
  output:
    tuple val(meta), val(meta.genome), path(bam_primary_out), emit: "primary_bam"
	tuple val(meta), val(meta.genome_spike), path(bam_spike_out), emit: "spike_bam"
	tuple val(meta), path(tsv_spike), emit: "spikes"
	path tsv_spike, emit: "report_align"

  script:
  bam_primary_out = "${rep_id}.bam"
  bam_spike_out = "${rep_id}_spike.bam"
  tsv_spike = "${rep_id}_spike.tsv"
  """
  bash ${params.dir_resources}/ngsutilsj bam-best \
    --stats ${tsv_spike} \
    --unsorted \
    ${bam_primary} \
    ${bam_spike} \
    -- \
    ${bam_primary_out} \
    ${bam_spike_out}
  """
}