# =============================================================================
# TriEval: random-sampling replications for the processed real data.
# Run from the repository root:
#   source("R/random_sample_replications.R")
# Prepares response/design subsets and starting values; does not run GVEM.
# =============================================================================
 
library(dplyr)
source("R/main.R")

# ARC is the benchmark supported by this repository's cleaning script.
# Add other names only when data/<benchmark>_resp.RData exists and uses
# the same context ordering as data/dummy_matrix.RData.
benchmarks <- c("arc")
sample_sizes <- c(200L, 600L, 1000L)
n_reps <- 50L
base_seed <- 2027
output_dir <- "results/random_sampling"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dummy_full <- as.matrix(readRDS("data/dummy_matrix.RData"))
N_full <- nrow(dummy_full)
stopifnot(!anyNA(dummy_full), all(is.finite(dummy_full)),
          all(sample_sizes > ncol(dummy_full)), all(sample_sizes <= N_full))

sampling_index <- list()

for (nm in benchmarks) {
  resp_full <- readRDS(file.path("data", paste0(nm, "_resp.RData")))
  item_cols <- grep("^Item[0-9]+$", names(resp_full), value = TRUE)
  stopifnot(length(item_cols) > 0L, "groupfull" %in% names(resp_full))
  resp_full <- resp_full[, c(item_cols, "groupfull"), drop = FALSE]
  Y_full <- as.matrix(resp_full[item_cols])
  stopifnot(!anyNA(resp_full$groupfull), !anyNA(Y_full),
            all(is.finite(Y_full)), all(Y_full %in% c(0, 1)))

  # As in the full-data estimator, within-group row k must correspond to
  # dummy_full row k. Do not filter or reorder groups independently.
  group_rows <- split(seq_len(nrow(resp_full)), resp_full$groupfull, drop = TRUE)
  stopifnot(all(lengths(group_rows) == N_full))

  for (n in sample_sizes) {
    folder <- file.path(output_dir, nm, paste0("random_", n))
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)

    for (rep_id in seq_len(n_reps)) {
      seed <- base_seed +
        ((match(nm, benchmarks) - 1L) * length(sample_sizes) +
           match(n, sample_sizes) - 1L) * n_reps + rep_id
      set.seed(seed)
      ids <- sort(sample.int(N_full, n, replace = FALSE))

      # Same context IDs for every model, retaining every item.
      selected_rows <- unlist(
        lapply(group_rows, function(rows) rows[ids]), use.names = FALSE
      )
      resp <- resp_full[selected_rows, , drop = FALSE]
      rownames(resp) <- NULL
      dummy_matrix <- dummy_full[ids, , drop = FALSE]
      group <- droplevels(as.factor(resp$groupfull))
      I <- nlevels(group)
      J <- length(item_cols)
      S <- nrow(dummy_matrix)

      # Stop rather than silently redraw a rank-deficient probability sample.
      stopifnot(all(table(group) == S), nrow(resp) == I * S,
                qr(cbind(1, dummy_matrix))$rank == ncol(dummy_matrix) + 1L)
      all_start <- start_value(resp, J = J)

      prepared_data <- list(
        resp = resp, dummy_matrix = dummy_matrix, group = group,
        I = I, J = J, S = S, N_full = N_full, all_start = all_start,
        context_ids = ids, source_rows = selected_rows,
        benchmark = nm, method = "random", replicate = rep_id, seed = seed
      )
      file <- file.path(folder, sprintf("rep_%03d_inputs.rds", rep_id))
      saveRDS(prepared_data, file)
      sampling_index[[length(sampling_index) + 1L]] <- data.frame(
        benchmark = nm, method = "random", size = S,
        replicate = rep_id, seed = seed, file = file
      )
    }
  }
}

write.csv(do.call(rbind, sampling_index),
          file.path(output_dir, "sampling_index.csv"), row.names = FALSE)
message("Saved ", length(sampling_index), " input bundles in ", output_dir)
