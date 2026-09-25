# =============================================================================
# requirements.R
# Installs all R packages required to run the code in this repository.
# Run this script once before running any other scripts.
#
# Usage:
#   Rscript requirements.R
# =============================================================================
 
required_packages <- c(
  "mirt",      # IRT model fitting (used in main.R)
  "dplyr",     # data manipulation
  "tidyr",     # data reshaping (pivot_wider)
  "parallel",  # parallel computation via mclapply (simulation.R)
  "ggplot2"    # plotting (simulation.R)
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  } else {
    message(paste0("'", pkg, "' is already installed."))
  }
}

invisible(lapply(required_packages, install_if_missing))

message("All required packages are installed.")
