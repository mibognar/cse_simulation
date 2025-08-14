# cse_simulation

Codebase to the paper:

Multiverse Simulation to Explore the Impact of Analytical Choices on Type I and Type II errors in a Reaction Time Study.

## Layout

```
congruentSeq/ # Local R package (Rcpp) to generate sequences
data/ # input/output data directoru (created/used by scripts)
sim_param_data.R # Simulation parameter grid + data generation
model.R # Model fitting code on simulated data
analysis.R # Extracting relevant data from the model objects and aggregating them
raw.R # (Optional) extracting raw CSE values to get a grasp of distribution
empirical_debrief.R # Extracting model parameters from the empirical datasets
batchtools.slurm.tmpl # SLURM template for future.batchtools
renv.lock # Locked dependency versions
.gitignore, session.txt (R session), etc.

```

## Prerequisites

- **R >= 4.2**
- **C++ toolchain** (for Rcpp)
- **Internet access for 'renv::restore()'**


## Environment

To reproduce computations, first clone the repository:

`git clone https://github.com/mibognar/cse_simulation -b hpc`

### Restore packages

```r

install.packages("renv")
renv::restore()

```

*note: Because of sudo restrictions in my environment, I bundled all of my local packages from the hpc environment to renv.lock, which is redundant.*

### Build and install the local Rcpp package congruentSeq

```r

install.packages("congruentSeq", repos = NULL, type = "source")

```

Then, you can load the package with the usual `library("congruentSeq")` call

## Runtime

To reproduce the workflow you should have SLURM in your environment, then

1. `sbatch sim_param_data.R`
2. (Wait to finish)
3. `sbatch model.R`
4. (Wait to finish)
5. `sbatch analysis.R`

*note: You should wait for each job to finish to start the next if you don't set up a slurm job dependency chain. This workflow needs a lot of space,
in our environment, the resulted directory size is over 9.2B.
The exact file outputs/paths are defined in the scripts: ensure data/ exists and that you have write permission.*

## Citation

If you use or adapt this code, please consider citing this repository and our paper.
