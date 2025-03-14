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

This pipeline is also is missing R packages `clugPac`, `ascFunctions`, and `ascCutNRun`. These are not necessary for any steps other than the final report, and they will be updated shortly.

## Execution

The Nextflow pipeline is executed within the `conda` environment using the command:
```
nextflow run cutNrun.nf \
  --sample_table <sample_table> \
  --dir_out <output/directory/path> \
  --control_epitope <background_epitope>
```

At minimum, three parameters should be provided to the nextflow pipeline, as well:

`-sample_table`: Filename of the input CSV file containing all sample details.

`-dir_out`: Output directory into which outputs will be published.

`-control_epitope`: Control epitope used for background signal, typically "IgG". This value must match at least one sample within the `sample_table` *per sample/condition combination*. All cell line/condition combinations are expected to have at least one background sample.

### Input table

The primary input of the `cutNrun.nf` pipeline is a CSV file containing one sample per row with annotations for each file and FastQ files for reads 1 and 2. The table format is:

project | name | cell_line | epitope | condition | replicate | R1 | R2 
--- | --- | --- | --- | --- | --- | --- |---
run1 | H358_MYC_WT_1 | H358      | MYC     | WT        |     1     | R1.fastq.gz | R2.fastq.gz
run1 | H358_MYC_WT_1 | H358      | MYC     | WT        |     2     | R1.fastq.gz | R2.fastq.gz
run1 | H358_IgG_WT_1 | H358      | IgG     | WT        |     1     | R1.fastq.gz | R2.fastq.gz

## Output

Output file structure is separated into `pooled` and `replicate` folders, with analyses run on individual replicates stored in the latter. Peak calling is performed by [`MACS3`](https://github.com/macs3-project/MACS). For individual replicates, no background is used: signal in a given region is compared to signal across the genome. Note that this is not ideal, particularly in samples with complex genomes. `MACS3` performs pooled peak calling by concatenating all replicates (treatment and background), then calling peaks using background samples to normalize for sequencing biases, copy number, etc. Note that this pipeline runs under the assumption that treatment and background samples should be produced within *the same sequencing run*, *the same cell type*, and *under the same conditions*. To modify this, alter the columns in the input CSV.

For instance, to use IgG background samples from one sequencing project as controls for another project, assign the replicates the same project name. This is not advised, however.
```
data
├── pooled
│   ├── H358_MYC_WG
│   │   ├── meme
│   │   │   ├── centrimo
│   │   │   ├── fimo
│   │   │   └── sea
│   │   ├── peaks
│   │   │   └── cuts
│   │   └── qc
└── replicates
    ├── H358_MYC_WT_1
    │   ├── align
    │   ├── peaks
    │   └── qc
    ├── H358_MYC_WT_2
    │   ├── align
    │   ├── peaks
    │   └── qc
    └── H358_IgG_WG_1
        ├── align
        ├── peaks
        └── qc
  
```
