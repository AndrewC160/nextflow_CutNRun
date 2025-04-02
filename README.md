# Cut&Run Nextflow pipeline

[Nextflow](https://github.com/nextflow-io/nextflow) is a popular scripting language that allows for tasks (defined as "processes") to be defined with inputs (dependencies) and outputs, and subsequent tasks that operate on these outputs will be run once all dependencies are ready. Tasks are parallellized, and a failed pipeline can be resumed following failure to pick up where it left off. In the event a script along the pipeline is changed or an input is altered, the process(es) that use that script/input will be repeated with the updated code/inputs. Nextflow is a widely used and readily available system that can be easily shared.

The `cutNrun.nf` script is the base script for the pipeline. When run, this script will subsequently call module scripts (stored in the `modules` directory) that contain specific steps for each process. Additional R scripts that will be called are included in the `R` folder, including a parameterized RMarkdown that will be used to generate an output document.

![*Example Cut&Run pipeline for KLF5 detection using an IgG background.*](https://github.com/user-attachments/assets/befeb93d-0720-43d1-aed0-74f950c133a8)

## Installation

The Cut&Run Nextflow pipeline is available on Github at [AndrewC160/nextflow_CutNRun](https://github.com/AndrewC160/nextflow_CutNRun). To install, clone this repository using:

`git clone https://github.com/AndrewC160/nextflow_CutNRun.git`

This repo includes a [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) environment stored in YML format, `conda_env.yml`. Presuming `conda` has been installed and is up-to-date, the `cutnrun` environment can be created using the command:

`conda env create -f conda_env.yml`

The environment can then be activated using:

`conda activate cutnrun`

**NOTE:** This version of the conda environment does not include [`ngsutilsj`](https://github.com/compgen-io/ngsutilsj), but a self-executing jar is provided within the pipeline environment. It should be added to the PATH once the `cutnrun` environment is activated. From within the `nextflow_cutNrun` directory, run:

`PATH=$(readlink -e ngsutilsj):$PATH`

This pipeline is also is missing R packages `clugPac`, `ascFunctions`, and `ascCutNRun`. These will be updated later, but for now essential functions are included in the R/R_functions directory which is sourced as necessary. All other R packages are included within the `conda_evn.yml` file.

## Configuration

By default Nextflow will operate within the bounds of the system it is run on (i.e. respecting memory/CPU limitations, etc.), but the nextflow.config file can be edited to set these limits explicitly.

### Resources

This pipeline requries several accessory files to run, and these are typically stored in a resources directory. This folder can be saved in the `nextflow_cutNrun` directory, or it can be specified using the `--dir_resources` argument. This folder should contain the following (though specific locations for each can be provided separately):

#### Bowtie2 indices

`--bt2_idx /path/to/bowtie2_index/hg38/hg38` (Primary genome, hg38 for instance)

`--bt2_spike /path/to/bowtie2_index/sac3/sac3` (Spike genome, Sac3 for instance)

Alignment is performed using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml), and so a Bowtie2 index is required for both the sample and spike genomes. [These can be generated manually](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome), but many pre-built indices are available for download ([here](https://benlangmead.github.io/aws-indexes/bowtie), for instance).

A Bowtie index should consist of a single directory *named after the genome*. Within this directory, Bowtie files should all have a prefix that matches the containing folder. For instance, the Bowtie2 index for the `hg38` genome is:
```
hg38
├── hg38.1.bt2
├── hg38.2.bt2
├── hg38.3.bt2
├── hg38.4.bt2
├── hg38.rev.1.bt2
└── hg38.rev.2.bt2
```

The index directories should be provided to the pipeline with the index prefix using the for the primary genome and for spike the genome.

#### Genome FASTA

`--fasta fasta/filename.fa`

Only required if MEME Suite functions are to be used.

#### Genome blacklists

--blacklist blacklist_file.bed

A bedfile of genomic blacklists should be provided. These are available for many genomes here, but keep in mind these files should be unzipped.

#### Seqsizes

--seqsizes hg38_seqsizes.tsv

TSV file with two columns: one for chromosome names, another for their length in basepairs. For instance:
```
chr1    248956422
chr2    242193529
chr3    198295559
chr4    190214555
```

#### Motif database

`--motif_db motif_database.meme`

Only required if [MEME Suite](https://meme-suite.org/meme/0) functions are used. MEME file of transcription factor binding motifs to search for, for instance [HOCOMOCO-v12](https://hocomoco12.autosome.org/downloads_v12).

#### Gene annotation GTF

`--gene_gtf hg38_genes.gtf`

GTF file of gene annotations, for instance from [UCSC GoldenPath](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz) This file should also be bgzipped and indexed.

## Execution

The cutNrun.nf script is the workflow for the pipeline. When run, this script will subsequently call module scripts (stored in the modules directory) that contain specific steps for each process. Additional R scripts that will be called are included in the R folder, including a parameterized RMarkdown that will be used to generate an output document.

The Nextflow pipeline is executed (within the active `conda` environment) using the command:
```
nextflow run cutNrun.nf \
  --sample_table <sample_table> \
  --dir_out <output/directory/path> \
  --control_epitope <background_epitope>
```

At minimum, three parameters should be provided to the nextflow pipeline, as well:

`--sample_table`: Filename of the input CSV file containing all sample details.

`--dir_out`: Output directory into which outputs will be published.

`--control_epitope`: Control epitope used for background signal, typically "IgG". This value must match at least one sample within the `sample_table` *per sample/condition combination*. All cell line/condition combinations are expected to have at least one background sample.

For the pooling process to run, each cell line/condition combination (pool) is expected to have at least one background sample; those without are only run through the replicate phase of the pipeline. It is also possible for background samples to apply to multiple pools: for instance, an IgG background in WT H358 cells can be used in one pool testing for KLF5 binding and another for H3K27 acetlyation.

### Sample table

The primary input of the `cutNrun.nf` pipeline is a CSV file containing one sample per row with annotations for each file and FastQ files for reads 1 and 2. The table format is:

project | name | cell_line | epitope | condition | replicate | R1 | R2 
--- | --- | --- | --- | --- | --- | --- |---
run1 | H358_MYC_WT_1 | H358      | MYC     | WT        |     1     | R1.fastq.gz | R2.fastq.gz
run1 | H358_MYC_WT_1 | H358      | MYC     | WT        |     2     | R1.fastq.gz | R2.fastq.gz
run1 | H358_IgG_WT_1 | H358      | IgG     | WT        |     1     | R1.fastq.gz | R2.fastq.gz

Each replicate should have a unique combination of project/cell_line/epitope/replicate IDs. Sample names are currently allowed, but output names for replicates and pools are constructed using `project/cell line/epitope/replicate/condition` tags:

> Replicates: `<cell_line>_<epitope>_<condition>_<replicate>`

> Pools: `<cell_line>_<epitope>_<condition>_<project>`

### ROSE (Optional)

Enable with `--run_rose true`.

By default, pooled `narrowPeak` data along with control/background BAM files are subjected to a custom abbreviated implementation of the [Rank Ordering of Super-enhancers](https://www.sciencedirect.com/science/article/pii/S221501612030385X) protocol designed to detect super-enhancers among signal peaks. Note that the usefulness of ROSE will depend on the epitope studied, and is typically applied to Mediator/H3K27Ac binding. ROSE results are stored in the pooled sample `ROSE` directory, and include a summary HTML, a TSV file of annotated super-enhancers, and an RDS file that contains a [SummarizedExperiment](https://www.bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) object that can be read into R directly:

```
library(SummarizedExperiment)               # Library for using SummarizedExperiments.
se           <- readRDS("rose_results.rds") # Load SummarizedExperiment.
raw_counts   <- assay(se,"counts")          # Access raw count matrix.
norm_counts  <- assay(se,"norm")            # Access normalized count matrix.
lib_info     <- colData(se)                 # Get library information (BAM files, read counts, weights,etc.).
region_info  <- rowData(se)                 # Get peak region information (Stitched regions, rank, scores, nearest genes, etc.).
se_cutoff    <- metadata(se)                # Super-enhancer area cutoff.
```

### MEME Suite (Optional)

Enable with `--run_meme true`.

Processes are included for three [MEME Suite](https://meme-suite.org/meme/) functions to detect DNA binding motifs among peak summits. Note that broadPeaks and often narrowPeaks can be large enough that many MEME functions are prohibitively slow. The functions currently implemented include:

-[SEA](https://meme-suite.org/meme/tools/sea): Simple Enrichment analysis; determines if any motifs are significantly enriched within peak summits relative to randomized sequences of the same FASTA.

-[CENTRIMO](https://meme-suite.org/meme/tools/centrimo): Determine if motifs have a significant preference for a location within peak summits.

-[FIMO](https://meme-suite.org/meme/tools/fimo): Finding Individual Motif Occurences; identifies DNA binding motifs within peak summits.

### Stub files (Optional)

To test/troubleshoot pipeline outputs, “stub” files can be generated: an optional process will truncate all raw FASTQ files to include a limited number of reads (limit specified via `--trunc_count`) that will process quickly and allow the entire pipeline to run to completion. For instance, to prepare a stub run using 1 million reads, provide the following two arguments:
```
--truncate_fastqs true \
--truncate_count 1000000
```

## Output

Output files are stored in the directory specified using `--dir_out`, within which `pooled` and `replicate` folders will be created. Analyses of pooled samples are stored in the former, individual replicates in the latter. Peak calling is performed by [MACS3](https://github.com/macs3-project/MACS). For individual replicates, no background is used: signal in a given region is compared to signal across the genome. This is not ideal, particularly in samples with complex genomes. For pooled peak calling, MACS3 concatenates all replicates (treatment and background), then calls peaks using background samples to normalize for sequencing biases, copy number, etc. Note that this pipeline runs under the assumption that *treatment and background samples were produced within the same sequencing run, the same cell type, and under the same conditions*. To modify this, alter the columns in the input CSV. For instance, to use IgG background samples from one sequencing project as controls for another, assign the replicates the same project name. This is not advised, however, as background samples should be produced alongside treatment samples.

Among the outputs for each pool, the `_spike.tsv` file includes read counts as well as the total number of reads aligned to the primary genome vs. the spike genome. This is useful when normalizing signal between samples. The `_report.html` file includes basic quality control metrics and plots of pileup among high-scoring peaks, and the `_file_summary.tsv` file contains output file locations.
```
data
├── pooled
│   ├── H358_MYC_WT_run1
│   │   ├── meme
│   │   │   ├── centrimo
│   │   │   ├── fimo
│   │   │   └── sea
│   │   ├── peaks
│   │   │   └── cuts
│   │   ├── qc
│   │   │   └── H358_MYC_WT_run1_report.html
│   │   ├── ROSE
│   │   |   ├── H358_MYC_WT_run1_ROSEoutput.rds
│   │   |   ├── H358_MYC_WT_run1_ROSEoutput.tsv
│   │   │   └── H358_MYC_WT_run1_ROSE.html
│   │   ├── H358_MYC_WT_spike.tsv
│   │   └── H358_MYC_WT_run1_file_summary.tsv
└── replicates
    ├── H358_MYC_WT_1
    │   ├── align
    │   ├── peaks
    │   └── qc
    ├── H358_MYC_WT_2
    │   ├── align
    │   ├── peaks
    │   └── qc
    └── H358_IgG_WT_1
        ├── align
        ├── peaks
        └── qc  
```

### nextflow_cutNrun/work

By default, Nextflow will generate a `work` directory in the working directory from which it is executed, typically the `nextflow_cutNrun` directory. These folders contains all subsequent steps in the pipeline, and to understand their application and the finer points of pipeline troubleshooting please see the [Nextflow documentation](https://www.nextflow.io/docs/latest/cache-and-resume.html). These files should not be necessary beyond troubleshooting and caching, and they can be safely deleted once a pipeline has been completed to free up disk space. Note that once removed (either via `rm -r` or `nextflow clean`), however, subsequent pipeline runs cannot be `-resume`d, and if run *the pipeline will repeat all steps and replace files in the `--out_dir` directory with fresh copies*. This is because the `work` directory files are tracked for caching, *not* the files in the output directory.
