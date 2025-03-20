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
params.bt2_idx = "${params.dir_bowtie}/hg38"
params.bt2_spike = "${params.dir_bowtie}/sac3"
params.fasta = "${params.dir_bowtie}/hg38/hg38.fa"
params.blacklist = "${params.dir_resources}/blacklists/hg38-blacklist.bed"
params.seqsizes = "${params.dir_bowtie}/hg38/hg38_seqsizes.tsv"
params.motif_db = "${params.dir_resources}/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"
params.gene_gtf = "${params.dir_resources}/Homo_sapiens.GRCh38.104.chr.tabix.gtf.gz"
params.gene_gtf_idx = "${params.gene_gtf}.tbi"

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
include { poolReport } from "${params.dir_modules}/mod_poolReport.nf"

workflow {
  // Read CSV.
  Channel.fromPath(params.sample_table)
    .splitCsv(header: true)
    .map { row -> [tuple(row.cell_line,row.cond,row.project).join("_"),
                   tuple(row.cell_line,row.epitope,row.cond,row.project).join("_"),
                   tuple(row.cell_line,row.epitope,row.cond).join("_"),
                   tuple(row.cell_line,row.epitope,row.cond,row.rep).join("_"),
                   row.epitope,file(row.R1),file(row.R2)] }
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
  alignment_hg38(trimming.out.trimmed,file(params.bt2_idx))
  
  // Alignment, sac3.
  alignment_sac3(trimming.out.trimmed,file(params.bt2_spike))
  
  // Genome assignment.
  ch_bams = alignment_hg38.out.aligned.join(alignment_sac3.out.aligned,by: 0..4)
  bamBest(ch_bams)
  
  // Combine replicate spike summaries.
  bamBest.out.spike
    .filter { !it[4].contains(params.control_epitope) }
    .map { row -> [row[1],row[2],row[5]] }
    .groupTuple(by:0..1)
    .map { id,condition,bam_files -> tuple(id,condition,bam_files.collect { file(it) }) }
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
    .filter { !it[4].contains(params.control_epitope) }
    .map { row -> row[0,1,3,5] }
    .join(by:0..2,
      peakCallingNarrow.out.narrowPeaks
        .map {row -> row[0,1,3,5] } )
    .set { ch_frip_reps }
  
  calculateFRiP(ch_frip_reps,"narrowPeaks")
  
  // Sample pooling.
  //  Split bams into treatment/background.
  //  Create an ID that comprises <project>_<cell_line>_<condition> and use it 
  //  to cross Treatment INTO background (required as *some backgrounds can be 
  //  used by multiple treatments*, and cross() is *not commutative*).
  readFiltering_hg38.out.filtered
    .map { row -> row[0..2,4,5] }
    .groupTuple(by:0..3)
    .set { ch_pooled_all }
  
  ch_pooled_all
    .filter { it[3].contains(params.control_epitope) }
    .map { row -> [row[0],row[4]] }
    .set { ch_pooled_ctrl }
  
  ch_pooled_all
    .filter { !it[3].contains(params.control_epitope) }
    .set { ch_pooled_test }
  
  ch_pooled_ctrl.cross(ch_pooled_test)
    .map { it -> it[1][0..2,4] + [it[0][1]] }
    .set { ch_pooled_bams }
  
  // Pooled peak calling.
  peakCallingNarrowPooled(ch_pooled_bams,file(params.blacklist),file(params.seqsizes))
  peakCallingBroadPooled(ch_pooled_bams,file(params.blacklist),file(params.seqsizes))
  
  // Cut points.
  if(params.run_cuts){
    peakCallingNarrowPooled.out.summits.join(peakCallingNarrowPooled.out.narrowPeaks,by:0..2)
      .map { row -> row[1,3,4]}
      .set { cut_regions }
    
    ch_pooled_bams
      .map{ row -> row[1,3] }
      .mix(
        ch_pooled_bams.map{ row -> row[1,4] }
      )
      .transpose()
      .set { cut_bams }
    
    cut_regions.cross(cut_bams)
      .map { it -> it[1][0,1] + it[0][1..2] }
      .set {ch_cuts}
      
    getCutPoints_summits(ch_cuts.map{ row -> tuple(row[0,1,2]) },params.dir_pool,"summit")
    getCutPoints_narrows(ch_cuts.map{ row -> tuple(row[0,1,3]) },params.dir_pool,"narrowPeak")
  }
  if(params.run_meme){
    // Retrieve peak sequences.
    getSequences_summits(peakCallingNarrowPooled.out.summits,params.fasta,"summits")
    getSequences_narrows(peakCallingNarrowPooled.out.narrowPeaks,params.fasta,"narrowPeaks")
    
    // SEA
    memeSEA(getSequences_summits.out.seqs,params.motif_db,"summits")
  
    
    // FIMO
    memeFIMO_summits(getSequences_summits.out.seqs,params.motif_db,"summits")
    //memeFIMO_narrows(getSequences_narrows.out.seqs,params.motif_db,"narrowPeaks")
  
    // CENTRIMO
    memeCENTRIMO(getSequences_summits.out.seqs,params.motif_db,"summits")
  }
  
  // Pool report.
  trimming.out.fastqc
    .map{ row -> row[1,5,6] }
    .groupTuple()
    .map{ sys_idx, fq1, fq2 -> tuple(sys_idx,fq1.collect { file(it) },fq2.collect { file(it) } ) }
    .set { rpt_fqc_trim }
  readFiltering_hg38.out.fastqc
    .map { row -> row[1,5] }
    .groupTuple()
    .map { id,fq_file -> tuple(id,fq_file.collect{ file(it) }) }
    .set { rpt_fqc_filt }
  readFiltering_hg38.out.fastqc
    .map { row -> row[1,5] }
    .groupTuple()
    .map { id, fqc_files -> tuple( id, fqc_files.collect { file(it) }) }
    .set { rpt_fqc_filt }
  spike_tsvs
    .map{ row -> row[0,2] }
    .map{ id,spike_files -> tuple(id,spike_files.collect { file(it) })}
    .set { rpt_spike }
  peakCallingNarrow.out.narrowPeaks
    .filter{ !it[4].contains(params.control_epitope) }
    .map { row -> row[1,5] }
    .groupTuple()
    .map { id,peak_files -> tuple(id,peak_files.collect { file(it) }) }
    .set{ rpt_npks_npks_rep }
  peakCallingNarrow.out.treatBDG
    .filter{ !it[4].contains(params.control_epitope) }
    .map { row -> row[1,5,6] }
    .groupTuple()
    .map { id,peak_files,indices -> tuple(id,peak_files.collect { file(it) }, indices.collect { file(it) }) }
    .set{ rpt_npks_bdgs_rep }
  peakCallingNarrow.out.ctrlBDG
    .filter{ !it[4].contains(params.control_epitope) }
    .map { row -> row[1,5,6] }
    .groupTuple()
    .map { id,peak_files,indices -> tuple(id,peak_files.collect { file(it) }, indices.collect { file(it) }) }
    .set{ rpt_npks_bdgs_ctrl_rep }
  peakCallingBroad.out.broadPeaks
    .filter{ !it[4].contains(params.control_epitope) }
    .map { row -> row[1,5] }
    .groupTuple()
    .map { id,peak_files -> tuple(id,peak_files.collect { file(it) }) }
    .set{ rpt_bpks_bpks_rep }
  calculateFRiP.out.frip
    .map { row -> row[1,3] }
    .groupTuple()
    .map { id,frip_files -> tuple(id,frip_files.collect { file(it) }) }
    .set{ rpt_frip }
  peakCallingNarrowPooled.out.narrowPeaks.map { row -> row[1,3] }.set{ rpt_npks_npks_pool }
  peakCallingNarrowPooled.out.summits.map { row -> row[1,3] }.set { rpt_npks_sums_pool }
  peakCallingNarrowPooled.out.treatBDG.map { row -> row[1,3,4]}.set { rpt_npks_bdgs_pool }
  peakCallingNarrowPooled.out.ctrlBDG.map { row -> row[1,3,4]}.set { rpt_npks_bdgs_ctrl_pool }
  peakCallingBroadPooled.out.broadPeaks.map { row -> row[1,3] }.set{ rpt_bpks_bpks_pool }
  
  rpt_fqc_trim
    .join(rpt_fqc_filt)
    .join(rpt_spike)
    .join(rpt_npks_npks_rep)
    .join(rpt_npks_bdgs_rep)
    .join(rpt_npks_bdgs_ctrl_rep)
    .join(rpt_bpks_bpks_rep)
    .join(rpt_frip)
    .join(rpt_npks_npks_pool)
    .join(rpt_bpks_bpks_pool)
    .join(rpt_npks_bdgs_pool)
    .join(rpt_npks_bdgs_ctrl_pool)
    .join(rpt_npks_sums_pool)
    .set { rpt_inputs }
  
  poolReport(
    file("${params.dir_R}/qc_pool.Rmd"),
    rpt_inputs,
    params.control_epitope,
    file(params.sample_table),
    params.dir_out,
    file("${params.dir_R}/R_functions/"),
    file(params.gene_gtf),
    file(params.gene_gtf_idx))
  
// ROSE

}