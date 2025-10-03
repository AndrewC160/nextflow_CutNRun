#!/usr/bin/env nextflow

/*
 * Cut&Run processing pipeline
 */
// Parameters.
// Inputs.
params.sample_table
params.dir_out
params.control_epitope = "IgG"
params.truncate_fastqs = false
params.truncate_count = 100000
params.run_rose = false
params.run_fimo = false
params.run_fimo_nPeaks = false
params.run_sea = false
params.run_sea_nPeaks = false
params.run_centrimo = false
params.run_memechip = false
params.pool_report = true

// Directories.
params.dir_modules = "${projectDir}/modules"
params.dir_R = "${projectDir}/R"
params.dir_resources = "${projectDir}/resources"
params.dir_bowtie = "${params.dir_resources}/bowtie2_indices"

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
include { alignment as alignment_primary } from "${params.dir_modules}/mod_alignment.nf"
include { alignment as alignment_spike } from "${params.dir_modules}/mod_alignment.nf"
include { bamBest } from "${params.dir_modules}/mod_bamBest.nf"
include { readFiltering as readFiltering_primary } from "${params.dir_modules}/mod_readFiltering.nf"
include { readFiltering as readFiltering_spike } from "${params.dir_modules}/mod_readFiltering.nf"
include { peakCallingNarrow } from "${params.dir_modules}/mod_peakCallingNarrow.nf"
include { peakCallingNarrowPooled } from "${params.dir_modules}/mod_peakCallingNarrowPooled.nf"
include { peakCallingBroad } from "${params.dir_modules}/mod_peakCallingBroad.nf"
include { peakCallingBroadPooled} from "${params.dir_modules}/mod_peakCallingBroadPooled.nf"
include { ROSE } from "${params.dir_modules}/mod_ROSE.nf"
include { calculateFRiP } from "${params.dir_modules}/mod_calculateFRiP.nf"
include { combineSpikes } from "${params.dir_modules}/mod_combineSpikes.nf"
include { getSequences as getSequences_summits } from "${params.dir_modules}/mod_getSequences.nf"
include { getSequences as getSequences_narrows } from "${params.dir_modules}/mod_getSequences.nf"
include { memeSEA } from "${params.dir_modules}/mod_memeSEA.nf"
include { memeFIMO as memeFIMO_summits} from "${params.dir_modules}/mod_memeFIMO.nf"
include { memeFIMO as memeFIMO_narrows} from "${params.dir_modules}/mod_memeFIMO.nf"
include { memeCENTRIMO } from "${params.dir_modules}/mod_memeCENTRIMO.nf"
include { memeChIP } from "${params.dir_modules}/mod_memeChIP.nf"
include { getCutPoints as getCutPoints_summits } from "${params.dir_modules}/mod_getCutPoints.nf"
include { getCutPoints as getCutPoints_narrows } from "${params.dir_modules}/mod_getCutPoints.nf"
include { poolReport } from "${params.dir_modules}/mod_poolReport.nf"

workflow {
  // Read CSV.
  Channel.fromPath(params.sample_table)
    .splitCsv(header: true)
	.map { row -> [row.project,row.name,row.cell_line,row.epitope,row.cond,row.rep,file(row.R1),file(row.R2)]
		replicate_name = "${row.cell_line}_${row.epitope}_${row.cond}_${row.rep}"
		pool_name = "${row.cell_line}_${row.epitope}_${row.cond}"
		rep_dir = "${params.dir_out}/${row.project}/replicates/${replicate_name}"
		pool_dir = "${params.dir_out}/${row.project}/${pool_name}"
		meta = [
			proj:row.project,
			name:row.name,
			cell_line:row.cell_line,
			epitope:row.epitope,
			condition:row.cond,
			rep:row.rep,
			replicate_name: replicate_name,
			pool_name: pool_name,
			rep_dir: rep_dir,
			pool_dir: pool_dir,
			genome: params.bt2_idx.substring(params.bt2_idx.lastIndexOf("/")+1),
			genome_spike: params.bt2_spike.substring(params.bt2_spike.lastIndexOf("/")+1),
			R1: row.R1,
			R2: row.R2
		]
		[meta,row.R1,row.R2]
	}
    .set { ch_input }
	
	
	// Truncate FastQs (OPTIONAL)
	// For use when running pipeline with stub dataset.
	(ch_trunc,ch_full) = params.truncate_fastqs
		? [ch_input,Channel.empty()]
		: [Channel.empty(),ch_input]

	truncateFastQs(ch_trunc,params.truncate_count)
		.mix(ch_full)
		.set { ch_fastqs }
	
	// Trimming.
	trimming(ch_fastqs)
	
	// Alignment, primary.
	alignment_primary(trimming.out.trimmed,file(params.bt2_idx))
		.aligned
		.map { meta,bam_file -> [meta.replicate_name,meta,bam_file]}
		.set{ ch_alg_primary }
	
	// Alignment, sac3.
	alignment_spike(trimming.out.trimmed,file(params.bt2_spike))
		.aligned
		.map { meta,bam_file -> [meta.replicate_name,meta,bam_file]}
		.set{ ch_alg_spike }
	
	// Genome assignment.
	ch_bams = ch_alg_primary.join(ch_alg_spike,by:0)

	bamBest(ch_bams)
	
	bamBest.out.spikes
		.map {meta,tsv_file -> [meta.epitope,meta.pool_name,meta.pool_dir,tsv_file]}
		.filter{ !it[0].contains(params.control_epitope) }
		.groupTuple(by:0..2)
		.set{ ch_pool_spikes }
	
	// Read filtering
	readFiltering_primary(bamBest.out.primary_bam,file(params.blacklist))
	readFiltering_spike(bamBest.out.spike_bam,file(params.blacklist))
	
	// Peak calling.
	peakCallingNarrow(readFiltering_primary.out.filtered,file(params.blacklist),file(params.seqsizes))
	peakCallingBroad(readFiltering_primary.out.filtered,file(params.blacklist),file(params.seqsizes))
	
	// FRiP scores.
	readFiltering_primary.out.filtered
		.map {meta,bam_file -> [meta.epitope,meta.replicate_name,bam_file]}
		.join(by:0..1,
			peakCallingNarrow.out.narrowPeaks
				.map { meta,bed_file -> [meta.epitope,meta.replicate_name,meta,bed_file] } )
		.set { ch_frip_reps }

	calculateFRiP(ch_frip_reps,"narrowPeaks")

	ch_input
		.map { meta,fq1,fq2 -> [
		meta.proj,
		meta.name,
		meta.cell_line,
		meta.epitope,
		meta.condition,
		meta.rep,
		meta.replicate_name,
		meta.pool_name,
		meta.rep_dir,
		meta.pool_dir,
		meta.genome,
		meta.genome_spike].join(",")}
		.collectFile(name:"reps.csv",newLine:true)
		.set { ch_reps_csv }

	// Sample pooling.
	// Split bams into treatment/background.
	//  Use cell line + condition to cross Treatment INTO Background (required 
	//  because backgrounds can be used by multiple treatments, and cross() 
	//  is not commutative).
	// NOTE: must sort groups within tuples; non-deterministic outputs means that filtering/
	//  grouping/re-combining creates different orders of files in each tuple that render
	//  subsequent process caches obsolete because their input tables differ.
	
	combineSpikes(ch_pool_spikes)
	
	readFiltering_primary.out.filtered
		.map{ meta,bam_file -> [
			 meta.cell_line,
			 meta.condition,
			 meta.proj,
			 meta.pool_name,
			 meta.epitope,
			 meta.pool_dir,
			 meta.genome,
			 meta.genome_spike,
			 meta.replicate_name,
			 bam_file]}
		.groupTuple(by:0..7)
		.set { ch_pooled_all }
	
	ch_pooled_all
		.filter { it[4].contains(params.control_epitope) }
		.map { row -> [row[0,1].join("_"),row[8],row[9]] }
		.set { ch_pooled_ctrl }
	
	ch_pooled_all
		.filter { !it[3].contains(params.control_epitope) }
		.map { row -> [row[0,1].join("_"),row] }
		.set { ch_pooled_test }
		
	ch_pooled_ctrl.cross(ch_pooled_test)
		.map { it -> it[1] + [it[0][1..2]] }
		.map { id,tr,ct -> [[
			proj:tr[2],
			pool_name:tr[3],
			cell_line:tr[0],
			epitope:tr[4],
			condition:tr[1],
			pool_dir:tr[5],
			genome:tr[6],
			genome_spike:tr[7],
			replicates:tr[8].sort().join(":"),
			replicates_bkg:ct[0].sort().join(":")],
			tr[9],ct[1]]}
		.set { ch_pooled_bams }

	// Pooled peak calling.
	peakCallingNarrowPooled(ch_pooled_bams,file(params.blacklist),file(params.seqsizes))
	peakCallingBroadPooled( ch_pooled_bams,file(params.blacklist),file(params.seqsizes))
	
	// Pooled peak FRiP scores.
	readFiltering_primary.out.filtered
		.map {meta,bam_file -> [meta.epitope,meta.replicate_name,bam_file]}
		.join(by:0..1,
			peakCallingNarrowPooled.out.narrowPeaks
				.map { meta,bed_file -> [meta.epitope,meta.replicate_name,meta,bed_file] } )
		.set { ch_frip_reps }
	
	// Retrieve peak sequences.
	getSequences_summits(peakCallingNarrowPooled.out.summits,params.fasta,"summits")
	getSequences_narrows(peakCallingNarrowPooled.out.narrowPeaks,params.fasta,"narrowPeaks")
	
	// MEME suite.
	//ch_fimo = Channel.empty()
	//ch_fimo_npks = Channel.empty()
	//ch_sea = Channel.empty()
	//ch_sea_npks = Channel.empty()
	//ch_centrimo = Channel.empty()	
	if(params.run_memechip){
		memeChIP(getSequences_summits.out.seqs,file(params.motif_db),"summits")
			.set {ch_memechip}
	}
	if(params.run_fimo){
		memeFIMO_summits(getSequences_summits.out.seqs,file(params.motif_db),"summits")
			.set {ch_fimo}
	}
	if(params.run_fimo_nPeaks){
		memeFIMO_narrows(getSequences_narrows.out.seqs,file(params.motif_db),"narrowPeaks")
			.set {ch_fimo_npks}
	}
	if(params.run_sea){
		memeSEA(getSequences_summits.out.seqs,file(params.motif_db),"summits")
			.set {ch_sea}
	}
	if(params.run_sea_nPeaks){
		memeSEA(getSequences_narrows.out.seqs,file(params.motif_db),"narrowPeaks")
			.set {ch_sea_npks}
	}
	if(params.run_centrimo){
		memeCENTRIMO(getSequences_summits.out.seqs,file(params.motif_db),"summits")
			.set {ch_centrimo}
	}

	// Rank Ordering of Super-Enhancers (ROSE).
	ch_rose = Channel.empty()
	if(params.run_rose){
		ROSE(peakCallingNarrowPooled.out.ROSE,
			file("${params.dir_R}/ROSE_asclab.Rmd"),
			file("${params.dir_R}/R_functions/"),
			file(params.gene_gtf),
			file(params.gene_gtf_idx))
			.set {ch_rose}
	}

	// Report.
	// Stage all reported files to directory to be detected/handled in R.
	ch_pooled_bams
		.map { meta,treat_bams,ctrl_bams -> [
		meta.proj,
		meta.pool_name,
		meta.cell_line,
		meta.epitope,
		meta.condition,
		meta.replicates,
		meta.replicates_bkg,
		meta.pool_dir,
		meta.genome,
		meta.genome_spike,
		treat_bams.name.join(":"),
		ctrl_bams.name.join(":")].join(",")}
		.collectFile(name:"pools.csv",newLine:true)
		.set { ch_pool_csv }
	
	poolReport(
		ch_input.map { meta,fq1,fq2 -> meta.proj }.first(),
		ch_reps_csv,
		ch_pool_csv,
		file(params.sample_table),
		params.control_epitope,
		params.dir_out,
		file("${params.dir_R}/R_functions/"),
		file(params.seqsizes),
		file(params.gene_gtf),
		file(params.gene_gtf_idx),
		file("${params.dir_R}/qc_report.Rmd"),
		trimming.out.report_qc.collect(),
		alignment_primary.out.report_qc.collect(),
		alignment_spike.out.report_qc.collect(),
		bamBest.out.report_align.collect(),
		readFiltering_primary.out.report_qc.collect(),
		readFiltering_primary.out.report_align.collect(),
		peakCallingNarrow.out.report_qc.collect(),
		peakCallingNarrow.out.report_peaks.collect(),
		peakCallingBroad.out.report_qc.collect(),
		peakCallingBroad.out.report_peaks.collect(),
		peakCallingNarrowPooled.out.report_qc.collect(),
		peakCallingNarrowPooled.out.report_peaks.collect(),
		peakCallingBroadPooled.out.report_qc.collect(),
		peakCallingBroadPooled.out.report_peaks.collect(),
		calculateFRiP.out.report_qc.collect(),
		combineSpikes.out.report_align.collect()
		)
}
