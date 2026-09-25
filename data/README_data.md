# Data

## Raw responses
 
The Benchmark responses (ARC Challenge as example here) are obtained from the 
public DOVE Lite dataset:

https://huggingface.co/datasets/nlphuji/DOVE_Lite

Save the raw response table as `data/output.csv`, then run
`source("R/real_data_clean.R")` from the repository root. The script expects
`sample_index`, `model`, `score`, and the five dimension columns named in its
column-selection step.

## Processed files

The preprocessing script writes:

| File | Contents |
|---|---|
| `arc_resp.RData` | Item responses and model-group labels |
| `dummy_matrix.RData` | Context design matrix |
| `arc_df_final_wide.RData` | Cleaned responses with context variables |

These are RDS files despite their `.RData` extensions; read them with `readRDS()`.
Processed responses are not included in this submission.

## Simulation coefficients

`main_alpha_results.csv` contains the model-specific prompt coefficients used
by `R/simulation.R`. Model columns use numeric labels. This table allows the
simulation study to run without downloading the raw responses.

New full-data fits write a coefficient table to
`results/real_data/main_alpha_results.csv`. Copy that file to
`data/main_alpha_results.csv` only if the simulation inputs are to be replaced.
