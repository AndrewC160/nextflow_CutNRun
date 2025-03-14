#!/usr/bin/env nextflow

/*
 * Parameters.
 */
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
    .map { row -> tuple(row.proj,row.name,row.cell_line,row.epitope,row.cond,row.rep,file(row.R1),file(row.R2)) }
    .set { ch_input }
  
  // Detect duplicate names.
  
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
  bt_idx = file(params.bt2_hg38)
  alignment_hg38(trimming.out.trimmed,bt_idx).aligned
  
  // Alignment, sac3.
  bt_idx = file(params.bt2_sac3)
  alignment_sac3(trimming.out.trimmed,bt_idx).aligned
  
  // Genome assignment.
  ch_bams = alignment_hg38.out.aligned.join(alignment_sac3.out.aligned,by: 0..5)
  bamBest(ch_bams)
  
  // Filtering.
  blk_lst = file(params.blacklist)
  readFiltering_hg38(bamBest.out.hg38,"hg38",blk_lst)
  readFiltering_sac3(bamBest.out.sac3,"sac3",blk_lst)
  
  // Peak calling.
  seq_szs = file(params.seqsizes)
  peakCallingNarrow(readFiltering_hg38.out.filtered,blk_lst,seq_szs,params.dir_reps)
  peakCallingBroad(readFiltering_hg38.out.filtered,blk_lst,seq_szs,params.dir_reps)
  
  // FRiP scores
  readFiltering_hg38.out.filtered
    .map { row -> tuple(row[1],row[6]) }
    .join(
      peakCallingNarrow.out.narrowPeaks 
        .map {row -> tuple(row[1],row[6]) } )
    .set { ch_frip_reps }
  calculateFRiP(ch_frip_reps,"narrowPeaks")
  
  // Sample pooling.
  //  Split bams into treatment/background.
  //  Create an ID that comprises <project>_<cell_line>_<condition> and use it 
  //  to cross Treatment INTO background (required as *some backgrounds can be 
  //  used by multiple treatments*, and cross() is *not commutative*).
  readFiltering_hg38.out.filtered
    .map { row -> tuple(row[0,2,4].join("_"),row[0],row[2..4].join("_"),row[2],row[3],row[4],row[6]) }
    .groupTuple(by:0..5)
    .set { ch_pooled_all }

  ch_pooled_all
    .filter { it[4].contains(params.control_epitope) }
    .map { row -> tuple(row[0],row[6]) }
    .set { ch_pooled_ctrl }
    
  ch_pooled_all
    .filter { !it[4].contains(params.control_epitope) }
    .set { ch_pooled_test }
  
  ch_pooled_ctrl.cross(ch_pooled_test)
    .map { it -> tuple(it[1][1],it[1][2],it[1][3],it[1][4],it[1][5],it[1][6],it[0][1]) }
    .set { ch_pooled_bams }
  
  // Pooled peak calling.
  peakCallingNarrowPooled(ch_pooled_bams,blk_lst,seq_szs)
  peakCallingBroadPooled(ch_pooled_bams,blk_lst,seq_szs)
  
  // Combine replicate spike summaries.
  bamBest.out.spike
    .filter { !it[3].contains(params.control_epitope) }
    .map { row -> tuple(row[2..4,0].join("_"),row[6]) }
    .groupTuple()
    .set { spike_tsvs }
  combineSpikes(spike_tsvs)
  
  // Retrieve peak sequences.
  getSequences_summits(peakCallingNarrowPooled.out.summits,params.fasta_hg38,"summits")
  getSequences_narrows(peakCallingNarrowPooled.out.narrowPeaks,params.fasta_hg38,"narrowPeaks")
  
  // Cut points.
  if(params.run_cuts){
    ch_pooled_bams
    .map{ row -> tuple(row[0],row[1],row[2],row[3],row[4],row[5]) }
    .mix(
      ch_pooled_bams.map{ row -> tuple(row[0],row[1],row[2],row[3],row[4],row[6]) }
    )
    .transpose() 
    .combine(  
      peakCallingNarrowPooled.out.summits.map{ row -> tuple(row[5]) }
    )
    .combine(  
      peakCallingNarrowPooled.out.narrowPeaks.map{ row -> tuple(row[5]) }
    )
    .set {ch_cuts}
    
    getCutPoints_summits(ch_cuts.map{ row -> tuple(row[1],row[5],row[6])},params.dir_pool,"summit")
    getCutPoints_narrows(ch_cuts.map{ row -> tuple(row[1],row[5],row[7])},params.dir_pool,"narrowPeak")
  }
  
  if(params.run_meme){
    // SEA
    memeSEA(getSequences_summits.out.seqs,params.motif_db,"summits")
    
    // FIMO
    memeFIMO_summits(getSequences_summits.out.seqs,params.motif_db,"summits")
    memeFIMO_narrows(getSequences_narrows.out.seqs,params.motif_db,"narrowPeaks")
    
    // CENTRIMO
    memeCENTRIMO(getSequences_summits.out.seqs,params.motif_db,"summits")
  }
  // ROSE
}