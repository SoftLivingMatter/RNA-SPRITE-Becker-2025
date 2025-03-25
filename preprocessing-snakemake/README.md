# Preprocessing Snakemake

> Part of RNA-SPRITE-Becker-2025 for determining RNA sequences clustered with SPRITE

## About
This workflow is adapted from the [sprite2.0-pipeline](https://github.com/GuttmanLab/sprite2.0-pipeline)
repository from the publication ["RNA promotes the formation of spatial compartments in the nucleus"](https://doi.org/10.1016/j.cell.2021.10.014).

Alterations include:
- Workflow structure follows best practices recommendations
- Added multiqc reporting for the entire dataset
- Added barcode template and simplified customization
- Included statistics of ligation efficiency and reporting
- Added fastqc and cluster deduplication
- Better slurm organization and environment usage
- Support for alignment with STAR or bowtie2

## Installation
Running the workflow requires [snakemake](https://snakemake.readthedocs.io/en/stable/index.html).
All work was performed with version 7.20.0.  All dependencies are constrained
in the workflow and requires a system installation of conda/mamba and singularity.
For slurm integration, you will need to use a profile or executor to work
with the scheduler.

## Usage
Configuration is set in two files.  The first, `paths.yaml` contains information
on all file locations and can be changed to suite individual needs.  Ensure
any modified file names contain the same wildcards.

The `parameters.yaml` file must be set prior to each run to match experimental
design and desired outputs/filters.
- All output files will be generated in the `workdir` directory.
- If set, `email` will attempt to email the provided address when an error occurs.
- The `input_fastq` is a dictionary to list multiple experiments to analyze.
Each entry needs a unique name and the location of the fastq.gz files to analyze
with wildcards for the `sample` and `read`.  The workflow will query the file system
to find samples and reads to analyze.
- The barcode ids can be provided as a text file, or if the file is not present,
a template and plate layout can be provided to generate the ids on demand.  We find
the plate layout is less error prone and can be modified in a spreadsheet application.
- `num_tags` sets the number of barcodes
- The first round of alignment uses bowtie2 and requires a reference with repeated
sequences.  The next round can be set as `star` or `bowtie` and requires the
corresponding reference location.  
- The remaining settings control parameters and options for various steps which
can be modified as needed.
- Finally, the last section has containers for setting the url and version of
several tools used throughout the workflow.


Once properly configured, you can run `snakemake` to generate the sequences.
For slurm integration, you will need to use a profile or executor to work
with the scheduler.  To use the specified environments, run with `--use-conda` and
`--use-singularity`.
