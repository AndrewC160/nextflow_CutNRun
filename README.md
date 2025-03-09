# Cut&Run Nextflow pipeline

[Nextflow](https://github.com/nextflow-io/nextflow) is a popular scripting language that allows for tasks (defined as "processes") to be defined with inputs (dependencies) and outputs, and subsequent tasks that operate on these outputs will be run once all dependencies are ready. Tasks are parallellized, and a failed pipeline can be resumed following failure to pick up where it left off. In the event a script along the pipeline is changed or an input is altered, the process(es) that use that script/input will be repeated with the updated code/inputs. Nextflow is a widely used and readily available system that can be easily shared.

The `cutNrun.nf` script is the base script for the pipeline. When run, this script will subsequently call module scripts (stored in the `modules` directory) that contain specific steps for each process. Additional R scripts that will be called are included in the `R` folder, including a parameterized RMarkdown that will be used to generate an output document.

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


