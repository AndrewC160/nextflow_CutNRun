# Cut&Run Nextflow pipeline

[Nextflow](https://github.com/nextflow-io/nextflow) is a popular scripting language that allows for tasks (defined as "processes") to be defined with inputs (dependencies) and outputs, and subsequent tasks that operate on these outputs will be run once all dependencies are ready. Tasks are parallellized, and a failed pipeline can be `-resume`d to pick up where it left off. In the event a script along the pipeline is changed or an input is altered, the process(es) that use that script/input will be repeated with the updated code/inputs. Nextflow is a widely used and readily available system that can be easily shared.

![cutNrun_nextflow_pipeline](https://github.com/user-attachments/assets/b610e71e-22db-477b-b462-14b069277dbf)

## Installation

The Cut&Run Nextflow pipeline is available on Github at [AndrewC160/nextflow_CutNRun](https://github.com/AndrewC160/nextflow_CutNRun). To install, clone this repository using:

`git clone https://github.com/AndrewC160/nextflow_CutNRun.git`

This repo includes a [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) environment stored in YML format, `conda_env.yml`. Presuming `conda` has been installed and is up-to-date, the `cutnrun` environment can be created using the command:

`conda env create -f conda_env.yml`

The environment can then be activated using:

`conda activate cutnrun`

**NOTE:** This version of the conda environment does not include [`ngsutilsj`](https://github.com/compgen-io/ngsutilsj), but a self-executing jar is provided within the pipeline environment. It should be added to the PATH once the `cutnrun` environment is activated. From within the `nextflow_cutNrun` directory, run:

`PATH=$(readlink -e ngsutilsj):$PATH`

This pipeline is also is missing R packages `clugPac`, `ascFunctions`, and `ascCutNRun`. These will be updated later, but for now essential functions are included in the R/R_functions directory.

## Configuration

By default Nextflow will operate within the bounds of the system it is run on (i.e. respecting memory/CPU limitations, etc.), but the `nextflow.config` file can be edited to set these limits explicitly. 

## Execution

The `cutNrun.nf` script is the workflow for the pipeline. When run, this script will subsequently call module scripts (stored in the `modules` directory) that contain specific steps for each process. Additional R scripts that will be called are included in the `R` folder, including a parameterized RMarkdown that will be used to generate an output document.

The Nextflow pipeline is executed within the `conda` environment using the command:
```
nextflow run cutNrun.nf \
  --sample_table <sample_table> \
  --dir_out <output/directory/path> \
  --control_epitope <background_epitope>
```

At minimum, three parameters should be provided to the nextflow pipeline:

`--sample_table`: Filename of the input CSV file containing all sample details.

`--dir_out`: Output directory into which outputs will be published.

`--control_epitope`: Control epitope used for background signal, typically "IgG". This value must match at least one sample within the `sample_table` *per sample/condition combination*. 

For the pooling process to run, each cell line/condition combination (pool) is expected to have at least one background sample; those without are only run through the replicate phase of the pipeline. It is also possible for background samples to apply to multiple pools: for instance, an IgG background in WT H358 cells can be used in one pool testing for KLF5 binding and another for H3K27 acetlyation. 

### Sample table

The primary input of the `cutNrun.nf` pipeline is a CSV file containing one sample per row with annotations for each file and FastQ files for reads 1 and 2. The table format is:

project | name | cell_line | epitope | condition | replicate | R1 | R2 
--- | --- | --- | --- | --- | --- | --- |---
run1 | H358_MYC_WT_1 | H358      | MYC     | WT        |     1     | R1.fastq.gz | R2.fastq.gz
run1 | H358_MYC_WT_2 | H358      | MYC     | WT        |     2     | R1.fastq.gz | R2.fastq.gz
run1 | H358_IgG_WT_1 | H358      | IgG     | WT        |     1     | R1.fastq.gz | R2.fastq.gz

Each replicate should have a unique combination of project/cell_line/epitope/replicate IDs. Sample names are currently allowed, but output names for pools are constructed using project/cell line/epitope/condition tags:

- Pool index: `<cell_line>_<epitope>_<condition>_<project>`

## Output

Output files are stored in the directory specified using `--dir_out`, within which `pooled` and `replicate` folders will be created. Analyses of pooled samples are stored in the former, individual replicates in the latter. Peak calling is performed by [`MACS3`](https://github.com/macs3-project/MACS). For individual replicates, no background is used: signal in a given region is compared to signal across the genome. This is not ideal, particularly in samples with complex genomes. For pooled peak calling, `MACS3` concatenates all replicates (treatment and background), then calls peaks using background samples to normalize for sequencing biases, copy number, etc. Note that this pipeline runs under the assumption that treatment and background samples should be produced within *the same sequencing run*, *the same cell type*, and *under the same conditions*. To modify this, alter the columns in the input CSV. For instance, to use IgG background samples from one sequencing project as controls for another project, assign the replicates the same project name. This is not advised, however, as background samples should be produced alongside treatment samples.

Among the outputs for each pool, the `_spike.tsv` file includes read counts as well as the total number of reads aligned to the primary genome vs. the spike genome. This is useful when normalizing signal between samples. The `_report.html` file includes basic quality control metrics as well as 
```
data
├── pooled
│   ├── H358_MYC_WT_run1
│   │   ├── meme
│   │   │   ├── centrimo
│   │   │   ├── fimo
│   │   │   └── sea
│   │   ├── peaks
│   │   │   └── cuts
│   │   ├── qc
│   │   ├── H358_MYC_WT_spike.tsv
│   │   └── H358_MYC_WT_run1_report.html
└── replicates
    ├── H358_MYC_WT_1
    │   ├── align
    │   ├── peaks
    │   └── qc
    ├── H358_MYC_WT_2
    │   ├── align
    │   ├── peaks
    │   └── qc
    └── H358_IgG_WT_1
        ├── align
        ├── peaks
        └── qc
```

### nextflow_cutNrun/work

By default, `Nextflow` will generate a `work` directory in the working directory from which it is executed, typically the `nextflow_cutNrun` directory. These folders contains all subsequent steps in the pipeline, and to understand the finer points of pipeline troubleshooting please see [the Nextflow documenation](https://www.nextflow.io/docs/latest/cache-and-resume.html). These files should not be necessary beyond troubleshooting and caching, and they can be deleted safely once a pipeline has been completed to free up disk space. Note that once removed (either via `rm -r` or `nextflow clean`), however, subsequent pipeline runs cannot be `-resume`d, and if run the pipeline will replace files in the `--out_dir` directory with fresh copies.
