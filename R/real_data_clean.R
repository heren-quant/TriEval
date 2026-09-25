# =============================================================================
# real_data_clean.R
# Cleans and reshapes the raw benchmark data into the format required
# by the GVEM model, and constructs the prompt factor dummy matrix.
#
# INPUT:
#   data/output.csv          — raw benchmark responses
#
# OUTPUT:
#   data/arc_df_final_wide.RData  — cleaned wide-format data frame
#   data/arc_resp.RData           — item response matrix with group labels
#   data/dummy_matrix.RData       — prompt factor dummy matrix (S x P)
# =============================================================================
 
rm(list = ls())

"%!in%" <- Negate("%in%")

library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# Load raw data
# The raw data is publicly available at:
#   https://huggingface.co/datasets/nlphuji/DOVE_Lite
# Download the file, rename it to output.csv, and place it in the data/ folder.
# -----------------------------------------------------------------------------
arc_df <- read.csv("data/output.csv")

# -----------------------------------------------------------------------------
# Step 1: Select and rename relevant columns
# -----------------------------------------------------------------------------
arc_df_clean <- arc_df %>%
  dplyr::select(sample_index, model, starts_with("dimension"), score) %>%
  dplyr::rename(
    D1 = dimensions_1..enumerator,
    D2 = dimensions_2..separator,
    D3 = dimensions_3..choices_order,
    D4 = dimensions_4..instruction_phrasing_text,
    D5 = dimensions_5..shots
  )

# -----------------------------------------------------------------------------
# Step 2: Filter out incomplete prompt conditions
# Only keep D1-D5 combinations with exactly 500 responses
# -----------------------------------------------------------------------------
arc_df_filtered <- arc_df_clean %>%
  mutate(group_id = cur_group_id(), .by = c(D1, D2, D3, D4, D5)) %>%
  group_by(group_id) %>%
  filter(n() == 500) %>%
  ungroup()

# -----------------------------------------------------------------------------
# Step 3: Sort and re-assign group IDs
# Sort by model then dimensions to ensure consistent ordering.
# Re-assign group_id so that for every model, the first D1-D5 combination
# is ID 1 (crucial for dummy matrix alignment).
# -----------------------------------------------------------------------------
arc_df_final <- arc_df_filtered %>%
  dplyr::arrange(model, D1, D2, D3, D4, D5) %>%
  dplyr::group_by(D1, D2, D3, D4, D5) %>%
  dplyr::mutate(group_id = cur_group_id()) %>%
  dplyr::ungroup()

# -----------------------------------------------------------------------------
# Step 4: Reshape to wide format (one row per model x prompt condition)
# -----------------------------------------------------------------------------
arc_df_final_wide <- arc_df_final %>%
  pivot_wider(names_from = sample_index, values_from = score)

colnames(arc_df_final_wide)[8:107] <- paste0("Item", 1:100)

# -----------------------------------------------------------------------------
# Step 5: Construct the prompt factor dummy matrix
# Extract unique D1-D5 combinations from the first model only,
# then generate the dummy matrix for these unique prompt conditions.
# -----------------------------------------------------------------------------
unique_dimensions <- arc_df_final_wide %>%
  dplyr::filter(model == unique(model)[1]) %>%
  dplyr::select(D1, D2, D3, D4, D5) %>%
  mutate(across(everything(), as.factor))

dummy_matrix_unique <- model.matrix(~ ., data = unique_dimensions)[, -1]

# -----------------------------------------------------------------------------
# Step 6: Construct item response data frame with group labels
# -----------------------------------------------------------------------------
arc_resp <- arc_df_final_wide %>%
  dplyr::select(starts_with("Item"), group_id, model) %>%
  mutate(groupfull = as.numeric(as.factor(model))) %>%
  dplyr::select(-model, -group_id)

# -----------------------------------------------------------------------------
# Save outputs
# -----------------------------------------------------------------------------
saveRDS(arc_df_final_wide,  file = "data/arc_df_final_wide.RData")
saveRDS(arc_resp,           file = "data/arc_resp.RData")
saveRDS(dummy_matrix_unique, file = "data/dummy_matrix.RData")