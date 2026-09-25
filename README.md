# TriEval 

R implementation of TriEval for estimating prompt effects and model- and
task-level variation in benchmark responses. Estimation uses regularized Gaussian
variational expectation–maximization (RGVEM).

## Files

| File | Purpose |
|---|---|
| `R/main.R` | RGVEM estimation functions |
| `R/real_data_clean.R` | Preprocess ARC responses and construct the design matrix |
| `R/real_data_run.R` | Fit the model to the processed ARC data |
| `R/random_sample_replications.R` | Prepare repeated random subsets of the processed data |
| `R/simulation.R` | Run the simulation study |
| `data/main_alpha_results.csv` | Coefficients used to generate simulation data |
| `data/README_data.md` | Data source and preprocessing instructions |
| `requirements.R` | Install package dependencies |

## Setup

Run commands from the repository root. Install dependencies with:

```r
source("requirements.R")
```

Dependencies are `mirt`, `dplyr`, `tidyr`, `parallel`, and `ggplot2`.
The preprocessing script requires dplyr 1.1.0 or later for `.by` syntax.

## Simulation

The included coefficient table supplies the simulation inputs. Raw benchmark
responses are not required for this step.

```r
source("R/simulation.R")
```

Outputs are written to `results/simulation/`.

## Real-data analysis

Download and prepare the raw data as described in `data/README_data.md`, then run:

```r
source("R/real_data_clean.R")
source("R/real_data_run.R")
```

The fitted coefficient table is written to
`results/real_data/main_alpha_results.csv`. The simulation script reads the
included table from `data/main_alpha_results.csv`; replace that input explicitly
if simulations are to use a newly fitted table.

## Random-sampling replications

After preprocessing, run:

```r
source("R/random_sample_replications.R")
```

The script prepares 50 random samples without replacement at each of three
sizes: 200, 600, and 1,000 contexts. Each replicate uses a distinct seed.
Selected context IDs are shared across all models and items within a replicate.
The default benchmark is ARC. Additional benchmarks require response files
with the same context ordering as `data/dummy_matrix.RData`.

Each sample is saved as
`results/random_sampling/<benchmark>/random_<size>/rep_<number>_inputs.rds`.
The file contains the response data, matching design matrix, dimensions,
starting values, selected IDs, and seed. `sampling_index.csv` lists the files.
Re-running the script with the same settings replaces the same outputs.
No estimation is performed by this script.

To load a sample for estimation:

```r
source("R/main.R")
selected <- readRDS("results/random_sampling/arc/random_600/rep_001_inputs.rds")
list2env(selected[c("resp", "dummy_matrix", "group", "I", "J", "S",
                   "N_full", "all_start")], envir = .GlobalEnv)
```

Use the RGVEM call in `R/real_data_run.R` with these inputs. Sourcing that entire
script instead reloads the full dataset. To preserve the relative penalty
strength of its full-data example, set `lambda = 30 * S / N_full` and
`c = 0.04 * S` in the sample fit.

## Indices

In the code, `i` indexes models, `j` indexes tasks, and `s` indexes contexts.
