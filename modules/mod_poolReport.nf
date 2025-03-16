#!/usr/bin/env nextflow

process poolReport {
  tag "${samp_name}"
  publishDir "${params.dir_reps}/${samp_name}_${proj}/qc", mode: 'copy', pattern: "*.html"
  
  input:
    tuple val(proj), val(samp_name), val(cell_line), val(epitope), val(cond), val(rep), path(bam_hg38_in), path(bam_sac3_in)
  
  output:
    path "*"
  
  script:
  bam_hg38_out = "${samp_name}_hg38.bam"
  bam_sac3_out = "${samp_name}_sac3.bam"
  tsv_spike = "${samp_name}_spike.tsv"
  """
  Rscript -e 'rmarkdown::render("${params.dir_R}/qc_pool.Rmd",output_format="html_document",
  """
  //pool_name: "OS526_K27Ac_NONE_Project_ACEH_EH005"
  //tsv_spike: "/N/p/asclab/ASC-cutNrun/25MAR09-tf_examples/data/pooled/OS526_K27Ac_NONE_Project_ACEH_EH005/OS526_K27Ac_NONE_Project_ACEH_EH005_spike.tsv"
  //bdg_pool: "/N/p/asclab/ASC-cutNrun/25MAR09-tf_examples/data/pooled/OS526_K27Ac_NONE_Project_ACEH_EH005/peaks/OS526_K27Ac_NONE_treat_pileup.bdg.gz"
  //bdgs_reps: "/N/p/asclab/ASC-cutNrun/25MAR09-tf_examples/data/replicates/OS526_NONE_K27Ac_1/peaks/OS526_NONE_K27Ac_1_treat_pileup.bdg.gz /N/p/asclab/ASC-cutNrun/25MAR09-tf_examples/data/replicates/OS526_NONE_K27Ac_2/peaks/OS526_NONE_K27Ac_2_treat_pileup.bdg.gz /N/p/asclab/ASC-cutNrun/25MAR09-tf_examples/data/replicates/OS526_NONE_K27Ac_3/peaks/OS526_NONE_K27Ac_3_treat_pileup.bdg.gz"
  //nPeaks_pool: "/N/p/asclab/ASC-cutNrun/25MAR09-tf_examples/data/pooled/OS526_K27Ac_NONE_Project_ACEH_EH005/peaks/OS526_K27Ac_NONE_peaks.narrowPeak"
  //bPeaks_pool: "/N/p/asclab/ASC-cutNrun/25MAR09-tf_examples/data/pooled/OS526_K27Ac_NONE_Project_ACEH_EH005/peaks/OS526_K27Ac_NONE_peaks.broadPeak"
  //summits_pool: "/N/p/asclab/ASC-cutNrun/25MAR09-tf_examples/data/pooled/OS526_K27Ac_NONE_Project_ACEH_EH005/peaks/OS526_K27Ac_NONE_summits.bed"
  //ctrl_epitope: "IgG"
}
