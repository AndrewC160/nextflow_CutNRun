#!/usr/bin/env nextflow

/*
 * Cut&Run processing pipeline
 */
// Parameters.
// Inputs.
params.sample_table
params.dir_out
params.control_epitope = "IgG"
params.truncate_fastqs = true
params.truncate_count = 100000
params.run_meme = false
params.run_cuts = false

// Directories.
params.dir_modules = "${projectDir}/modules"
params.dir_R = "${projectDir}/R"
params.dir_resources = "${projectDir}/resources"
params.dir_bowtie = "${params.dir_resources}/bowtie2_indices"
params.dir_reps = "${params.dir_out}/replicates"
params.dir_pool = "${params.dir_out}/pooled"

// Accessory files.
params.bt2_hg38 = "${params.dir_bowtie}/hg38"
params.bt2_sac3 = "${params.dir_bowtie}/sac3"
params.fasta_hg38 = "${params.dir_bowtie}/hg38/hg38.fa"
params.blacklist = "${params.dir_resources}/blacklists/hg38-blacklist.bed"
params.seqsizes = "${params.dir_bowtie}/hg38/hg38_seqsizes.tsv"
params.motif_db = "${params.dir_resources}/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"

// Modules.
include { truncateFastQs } from "${params.dir_modules}/mod_truncateFastQs.nf"
include { trimming } from "${params.dir_modules}/mod_trimming.nf"
include { alignment as alignment_hg38 } from "${params.dir_modules}/mod_alignment.nf"
include { alignment as alignment_sac3 } from "${params.dir_modules}/mod_alignment.nf"
include { bamBest } from "${params.dir_modules}/mod_bamBest.nf"
include { readFiltering as readFiltering_hg38 } from "${params.dir_modules}/mod_readFiltering.nf"
include { readFiltering as readFiltering_sac3 } from "${params.dir_modules}/mod_readFiltering.nf"
include { peakCallingNarrow } from "${params.dir_modules}/mod_peakCallingNarrow.nf"
include { peakCallingNarrowPooled } from "${params.dir_modules}/mod_peakCallingNarrowPooled.nf"
include { peakCallingBroad } from "${params.dir_modules}/mod_peakCallingBroad.nf"
include { peakCallingBroadPooled} from "${params.dir_modules}/mod_peakCallingBroadPooled.nf"
include { calculateFRiP } from "${params.dir_modules}/mod_calculateFRiP.nf"
include { combineSpikes } from "${params.dir_modules}/mod_combineSpikes.nf"
include { getSequences as getSequences_summits } from "${params.dir_modules}/mod_getSequences.nf"
include { getSequences as getSequences_narrows } from "${params.dir_modules}/mod_getSequences.nf"
include { memeSEA } from "${params.dir_modules}/mod_memeSEA.nf"
include { memeFIMO as memeFIMO_summits} from "${params.dir_modules}/mod_memeFIMO.nf"
include { memeFIMO as memeFIMO_narrows} from "${params.dir_modules}/mod_memeFIMO.nf"
include { memeCENTRIMO } from "${params.dir_modules}/mod_memeCENTRIMO.nf"
include { getCutPoints as getCutPoints_summits } from "${params.dir_modules}/mod_getCutPoints.nf"
include { getCutPoints as getCutPoints_narrows } from "${params.dir_modules}/mod_getCutPoints.nf"

workflow {
  // Read CSV.
  Channel.fromPath(params.sample_table)
    .splitCsv(header: true)
    .map { row -> [tuple(row.cell_line,row.cond,row.project).join("_"),
                   tuple(row.cell_line,row.epitope,row.cond,row.project).join("_"),
                   tuple(row.cell_line,row.epitope,row.cond).join("_"),
                   row.project,row.name,row.cell_line,row.epitope,row.cond,row.rep,file(row.R1),file(row.R2)] }
    .set { ch_input }
  
  // Indices are used for grouping explicitly; so far they include:
  // System index: <cell_line>_<condition>_<project>
  // Sample index: <cell_line>_<epitope>_<condition>_<project>
  // Sample name: <cell_line>_<epitope>_<condition>
  
  // Truncate FastQs (OPTIONAL)
  // For use when running pipeline with stub dataset.
  (ch_trunc,ch_full) = params.truncate_fastqs
    ? [ch_input,Channel.empty()]
    : [Channel.empty(), ch_input]
  
  truncateFastQs(ch_trunc,params.truncate_count)
    .mix(ch_full)
    .set { ch_fastqs }
  
  // Trimming.
  trimming(ch_fastqs)
  
  // Alignment, hg38.
  alignment_hg38(trimming.out.trimmed,file(params.bt2_hg38))
  
  // Alignment, sac3.
  alignment_sac3(trimming.out.trimmed,file(params.bt2_sac3))
  
  // Genome assignment.
  ch_bams = alignment_hg38.out.aligned.join(alignment_sac3.out.aligned,by: 0..8)
  bamBest(ch_bams)
  
  // Combine replicate spike summaries.
  bamBest.out.spike
    .filter { !it[6].contains(params.control_epitope) }
    .map { row -> [row[1],row[2],row[9]] }
    .groupTuple(by:0..1)
    .set { spike_tsvs }

  combineSpikes(spike_tsvs)
  
  // Filtering.
  readFiltering_hg38(bamBest.out.hg38,"hg38",file(params.blacklist))
  readFiltering_sac3(bamBest.out.sac3,"sac3",file(params.blacklist))
  
  // Peak calling.
  peakCallingNarrow(readFiltering_hg38.out.filtered,file(params.blacklist),file(params.seqsizes),params.dir_reps)
  peakCallingBroad(readFiltering_hg38.out.filtered,file(params.blacklist),file(params.seqsizes),params.dir_reps)
  
  // FRiP scores
  readFiltering_hg38.out.filtered
    .filter { !it[6].contains(params.control_epitope) }
    .map { row -> row[1,4,9] }
    .join(by:0..1,
      peakCallingNarrow.out.narrowPeaks
        .map {row -> row[1,4,9] } )
    .set { ch_frip_reps }
  calculateFRiP(ch_frip_reps,"narrowPeaks")
  
  // Sample pooling.
  //  Split bams into treatment/background.
  //  Create an ID that comprises <project>_<cell_line>_<condition> and use it 
  //  to cross Treatment INTO background (required as *some backgrounds can be 
  //  used by multiple treatments*, and cross() is *not commutative*).
  readFiltering_hg38.out.filtered
    .map { row -> row[0..3,5..7,9] }
    .groupTuple(by:0..6)
    .set { ch_pooled_all }
  
  ch_pooled_all
    .filter { it[5].contains(params.control_epitope) }
    .map { row -> [row[0],row[7]] }
    .set { ch_pooled_ctrl }
  
  ch_pooled_all
    .filter { !it[5].contains(params.control_epitope) }
    .set { ch_pooled_test }
  
  ch_pooled_ctrl.cross(ch_pooled_test)
    .map { it -> it[1][0..7] + [it[0][1]] }
    .set { ch_pooled_bams }
  
  // Pooled peak calling.
  peakCallingNarrowPooled(ch_pooled_bams,file(params.blacklist),file(params.seqsizes))
  peakCallingBroadPooled(ch_pooled_bams,file(params.blacklist),file(params.seqsizes))
  
  // Cut points.
  if(params.run_cuts){
    ch_pooled_bams
      .map{ row -> row[0..7] }
      .mix(
        ch_pooled_bams.map{ row -> row[0..6,8] }
      )
      .transpose() 
      .combine(  
        peakCallingNarrowPooled.out.summits.map{ row -> tuple(row[7]) }
      )
      .combine(  
        peakCallingNarrowPooled.out.narrowPeaks.map{ row -> tuple(row[7]) }
      )
      .set {ch_cuts}
    getCutPoints_summits(ch_cuts.map{ row -> tuple(row[1],row[7],row[8]) },params.dir_pool,"summit")
    getCutPoints_narrows(ch_cuts.map{ row -> tuple(row[1],row[7],row[9]) },params.dir_pool,"narrowPeak")
  }
  
  if(params.run_meme){
    // Retrieve peak sequences.
    getSequences_summits(peakCallingNarrowPooled.out.summits,params.fasta_hg38,"summits")
    getSequences_narrows(peakCallingNarrowPooled.out.narrowPeaks,params.fasta_hg38,"narrowPeaks")
    
    // SEA
    memeSEA(getSequences_summits.out.seqs,params.motif_db,"summits")
    
    // FIMO
    memeFIMO_summits(getSequences_summits.out.seqs,params.motif_db,"summits")
    //memeFIMO_narrows(getSequences_narrows.out.seqs,params.motif_db,"narrowPeaks")
    
    // CENTRIMO
    memeCENTRIMO(getSequences_summits.out.seqs,params.motif_db,"summits")
  }
  // Pool report.
  trimming.out.fastqc.map { row -> row[1,9,10] }.set { rpt_trim }
  readFiltering_hg38.out.fastqc.map { row -> row[1,9] }.set { rpt_fqc_filt }
  spike_tsvs.map{ row -> row[0,2] }.set { rpt_spike }
  peakCallingNarrow.out.narrowPeaks.filter{ !it[6].contains(params.control_epitope) }.map { row -> row[1,9] }.groupTuple().set{ rpt_npks_npks_rep }
  peakCallingNarrow.out.treatBDG.map { row -> row[1,9,10]}.set { rpt_npks_bdgs_rep }
  peakCallingBroad.out.broadPeaks.filter{ !it[6].contains(params.control_epitope) }.map { row -> row[1,9] }.groupTuple().set{ rpt_bpks_bpks_rep }
  calculateFRiP.out.frip.map { row -> row[0,2] }.set{ rpt_frip }
  peakCallingNarrowPooled.out.narrowPeaks.map { row -> row[1,7] }.set{ rpt_npks_npks_pool }
  peakCallingNarrowPooled.out.summits.map { row -> row[1,7] }.set { rpt_npks_sums_pool }
  peakCallingNarrowPooled.out.treatBDG.map { row -> row[1,7,8]}.set { rpt_npks_bdgs_pool }
  peakCallingBroadPooled.out.broadPeaks.map { row -> row[1,7] }.set{ rpt_bpks_bpks_pool }
  rpt_trim
    .join(rpt_fqc_filt)
    .join(rpt_spike)
    .join(rpt_npks_npks_rep)
    .join(rpt_npks_bdgs_rep)
    .join(rpt_bpks_bpks_rep)
    .join(rpt_frip)
    .join(rpt_npks_npks_pool)
    .join(rpt_bpks_bpks_pool)
    .join(rpt_npks_bdgs_pool)
    .join(rpt_npks_sums_pool)
    .view()
  //trimming.out.fastqc.map { row -> row[1,9,10] }.view()
  // [OS526_IgG_RUNX2KO_ucdavis_23MAR06, /N/p/asclab/ASC-cutNrun/nextflow/nextflow_CutNRun/work/4e/9581562a8c0ad7d133dfcb6ea554b4/OS526_IgG_RUNX2KO_1_val_1_fastqc.zip, /N/p/asclab/ASC-cutNrun/nextflow/nextflow_CutNRun/work/4e/9581562a8c0ad7d133dfcb6ea554b4/OS526_IgG_RUNX2KO_1_val_2_fastqc.zip]
  //readFiltering_hg38.out.fastqc.map { row -> row[1,9] }.view()
  // [OS526_IgG_RUNX2KO_ucdavis_23MAR06, /N/p/asclab/ASC-cutNrun/nextflow/nextflow_CutNRun/work/1a/bb4537e81b27a3b3b680ad27abda7b/OS526_IgG_RUNX2KO_1_hg38_short_fastqc.zip]
  //spike_tsvs.view()
  // [OS526_K27Ac_NONE_CnR_H358_ELF, OS526_K27Ac_NONE, [/N/p/asclab/ASC-cutNrun/nextflow/nextflow_CutNRun/work/8f/82b8a1e24db69c769a29c482c84bdd/OS526_K27Ac_NONE_1_spike.tsv]]
  //peakCallingNarrow.out.narrowPeaks.filter{ !it[6].contains(params.control_epitope) }.map { row -> row[1,9] }.groupTuple().view()
  // [OS526_K27Ac_KO_Project_ACEH_EH005, [/N/p/asclab/ASC-cutNrun/nextflow/nextflow_CutNRun/work/bf/7b64e3f92ec3319c41ae22b7dea10a/OS526_K27Ac_KO_1_peaks.narrowPeak, /N/p/asclab/ASC-cutNrun/nextflow/nextflow_CutNRun/work/9a/0ba429f82fa4785df4a557a3cdb4fb/OS526_K27Ac_KO_2_peaks.narrowPeak]]
  //peakCallingBroad.out.broadPeaks.filter{ !it[6].contains(params.control_epitope) }.map { row -> tuple(row[1,9]) }.groupTuple().view()
  // [OS526_K27Ac_NEG_ucdavis_23MAR06, [/N/p/asclab/ASC-cutNrun/nextflow/nextflow_CutNRun/work/79/53a84498b9958cca3f8500ff46d398/OS526_K27Ac_NEG_2_peaks.broadPeak, /N/p/asclab/ASC-cutNrun/nextflow/nextflow_CutNRun/work/ff/7d14131d877e942b1c02b7903344a0/OS526_K27Ac_NEG_1_peaks.broadPeak]]
  //calculateFRiP.out.frip.map { row -> row[0,2] }.view()
  // [OS526_RUNX2KO_ucdavis_23MAR06, /N/p/asclab/ASC-cutNrun/nextflow/nextflow_CutNRun/work/99/d5919ddde2aff4e957aedef40bf2af/OS526_K27Ac_RUNX2KO_1_narrowPeaks_FRiP.txt]
  //peakCallingNarrowPooled.out.narrowPeaks.map { row -> row[1,7] }.view()
  // [OS526_K27Ac_NONE_CnR_H358_ELF, /N/p/asclab/ASC-cutNrun/nextflow/nextflow_CutNRun/work/cb/2bc2a8612c6a1bb3eb08087e6ebc97/OS526_K27Ac_NONE_peaks.narrowPeak]
  //peakCallingNarrowPooled.out.treatBDG.view()
  // [OS526_KO_Project_ACEH_EH005, OS526_K27Ac_KO_Project_ACEH_EH005, OS526_K27Ac_KO, Project_ACEH_EH005, OS526, K27Ac, KO, /N/p/asclab/ASC-cutNrun/nextflow/nextflow_CutNRun/work/94/c3612cd9ac95743ba20b762e654e0f/OS526_K27Ac_KO_treat_pileup.bdg.gz, /N/p/asclab/ASC-cutNrun/nextflow/nextflow_CutNRun/work/94/c3612cd9ac95743ba20b762e654e0f/OS526_K27Ac_KO_treat_pileup.bdg.gz.tbi]
  //peakCallingNarrowPooled.out.summits.map { row -> row[1,7] }.view()
  // [OS526_K27Ac_NONE_CnR_H358_ELF, /N/p/asclab/ASC-cutNrun/nextflow/nextflow_CutNRun/work/cb/2bc2a8612c6a1bb3eb08087e6ebc97/OS526_K27Ac_NONE_summits.bed]
  //peakCallingBroadPooled.out.broadPeaks.map { row -> row[1,7] }.view()
  // [OS526_K27Ac_RUNX2KO_ucdavis_23MAR06, /N/p/asclab/ASC-cutNrun/nextflow/nextflow_CutNRun/work/e0/a08ffcc516b50d66badbd196dac6d5/OS526_K27Ac_RUNX2KO_peaks.broadPeak]
  //pool_name: "OS526_K27Ac_NONE_Project_ACEH_EH005"
  //tsv_spike: "/N/p/asclab/ASC-cutNrun/25MAR09-tf_examples/data/pooled/OS526_K27Ac_NONE_Project_ACEH_EH005/OS526_K27Ac_NONE_Project_ACEH_EH005_spike.tsv"
  //bdg_pool: "/N/p/asclab/ASC-cutNrun/25MAR09-tf_examples/data/pooled/OS526_K27Ac_NONE_Project_ACEH_EH005/peaks/OS526_K27Ac_NONE_treat_pileup.bdg.gz"
  //bdgs_reps: "/N/p/asclab/ASC-cutNrun/25MAR09-tf_examples/data/replicates/OS526_NONE_K27Ac_1/peaks/OS526_NONE_K27Ac_1_treat_pileup.bdg.gz /N/p/asclab/ASC-cutNrun/25MAR09-tf_examples/data/replicates/OS526_NONE_K27Ac_2/peaks/OS526_NONE_K27Ac_2_treat_pileup.bdg.gz /N/p/asclab/ASC-cutNrun/25MAR09-tf_examples/data/replicates/OS526_NONE_K27Ac_3/peaks/OS526_NONE_K27Ac_3_treat_pileup.bdg.gz"
  //nPeaks_pool: "/N/p/asclab/ASC-cutNrun/25MAR09-tf_examples/data/pooled/OS526_K27Ac_NONE_Project_ACEH_EH005/peaks/OS526_K27Ac_NONE_peaks.narrowPeak"
  //bPeaks_pool: "/N/p/asclab/ASC-cutNrun/25MAR09-tf_examples/data/pooled/OS526_K27Ac_NONE_Project_ACEH_EH005/peaks/OS526_K27Ac_NONE_peaks.broadPeak"
  //summits_pool: "/N/p/asclab/ASC-cutNrun/25MAR09-tf_examples/data/pooled/OS526_K27Ac_NONE_Project_ACEH_EH005/peaks/OS526_K27Ac_NONE_summits.bed"
  //ctrl_epitope: "IgG"
// ROSE

}