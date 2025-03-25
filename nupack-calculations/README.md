# Nupack Calculations

> Part of RNA-SPRITE-Becker-2025 for calculating RNA free energy

## About
The nupack calculations workflow is a thin wrapper around two python scripts
for calculating the free energy of folding for monomers and dimers.  It aims
to optimize wall clock time of slurm cluster usage by generating batches of
sequences and requesting resources based on prior calls.  Initial profiling
showed that most batches required time and memory on the order of `O(n^2)` with
the longest sequence in a batch of monomers or total sequence length for monomers.

The `workflow` and `config` folders are used by `snakemake` to run new free
energy estimates.  The scripts and notebooks for optimizing performance on a
different cluster are provided in the `profiling` directory.

## Workflow
### Installation
Running the workflow requires [snakemake](https://snakemake.readthedocs.io/en/stable/index.html).
All work was performed with version 7.20.0.  To execute the free energy estimate,
an environment with nupack (v4.0.1.8) installed is also required.  Due to the licensing,
the user is required to create this environment themselves and include biopython (v1.78).
Detailed instructions for installing nupack can be found [on their website](https://docs.nupack.org/start/#installation-requirements).

### Usage
All configuration is set in the `config/config.yaml` file:

- All output files will be generated in the `workdir` directory.
- Set the location of the conda environment with `nupack` in the `nupack_env` variable.
- Multiple experiments can be listed under the `experiments:` dictionary.
At a minimum, an experiment needs a unique name, one estimate, and an input fasta.
For `monomer` estimates, free energy is calculated for all sequences in the `input_fasta`.
Optionally, setting the `monomer_csv` will limit calculations to those sequences
in the csv.
With `dimer` estimates, free energy is determined for any dimers specified in
the `paired_csv`.  Sequences in the `paired_csv` must be present in the `input_fasta.`
- The paths are set to reasonable defaults, but if you desire a different output
file tree you can update the strings as long as the same wildcards are present.

Once properly configured, you can run `snakemake` to generate the sequences.
For slurm integration, you will need to use a profile or executor to work
with the scheduler.  The resource estimates are set in the `workflow/Snakefile`
and can be updated in the `pred_time` and `pred_rss` functions.  To enable
progressive increases in resource usage, run snakemake with `--retries 3`.

The scripts used for calculating free energy can be found in the `workflow/scripts`
directory.

## Profiling
If you need to optimize resource usage for your own hardware and sequences, you
must first run a workflow to completion.  One option would be to increase the resource
estimates by a factor of say 10, and change the dependence on attempts to scale exponentially.
Running with `--retries 3` should eventually cause all or most jobs to complete.

### Installation
The provided scripts and notebook require [reportseff](https://github.com/troycomi/reportseff/)
for reading slurm usage, along with jupyter, pandas, numpy, and seaborn.

### Usage
Update the base directory for `profiling/get_usage.sh` to the slurm logs of
the test workflow.  This script will use reportseff to read all slurm usage
of successfully completed jobs and pass that output to a script to read the
lengths of each sequence.  The final report is used by the notebook in 
`profiling/analysis.ipynb` to generate fitted parameters.

The analysis notebook shows rendered results from our optimizations, indicating that
batching sequences into random groups of 10 provided the most consistent results.
Fitting to a quadratic model and scaling the output efficiency to 80% provided
a good balance between failed jobs and efficient usage.

As a final optimization, jobs that require more memory request more cores and scale
runtimes accordingly.  This keeps a constant fraction of a node occupied while
providing an asymptotic runtime of around 3 hours.  Sequences longer than 
40,000 are discarded as they consume too many resources.
