# =========================================================================
# Summarize and Aggregate SCM Simulation Results
# Last updated: 2025-10-14
# =========================================================================

# Load required packages
library(dplyr)
library(readr)

# Set the primary working directory for the project.
setwd("D:/R_work/Inflorescence")

# -------------------------------------------------------------------------
# Part 1: Parameters and Setup
# -------------------------------------------------------------------------
# Directory containing the individual SCM simulation results
scm_dir <- "Evolution_Model/SCM_Simulations_HR2_ARD"
# Path for the final aggregated output file
output_file <- "Evolution_Model/SCM_Summary_HR2_ARD.csv"

# List all CSV files that match the simulation output format
scm_files <- list.files(scm_dir, pattern = "^scm_summary_seed_\\d{4}\\.csv$", full.names = TRUE)
n_files <- length(scm_files)

if (n_files == 0) {
  stop("No SCM simulation files found in the specified directory: ", scm_dir)
}

# -------------------------------------------------------------------------
# Part 2: Read and Combine All Simulation Files
# -------------------------------------------------------------------------
scm_data_list <- list()

cat("Starting to read", n_files, "simulation result files...\n")

for(i in 1:n_files) {
  file_path <- scm_files[i]
  # Extract the 4-digit seed number from the filename
  seed_id <- gsub(".*seed_(\\d{4})\\.csv$", "\\1", basename(file_path))
  
  tryCatch({
    # Read the CSV and add the seed_id as a new column
    data <- read_csv(file_path, col_types = cols()) # Use read_csv for speed
    data$seed_id <- seed_id
    scm_data_list[[i]] <- data
    
    # Print progress periodically
    if(i %% 100 == 0 || i == n_files) {
      cat("Read", i, "/", n_files, "files\n")
    }
  }, error = function(e) {
    cat("Warning: Failed to read file", basename(file_path), ". Error:", e$message, "\n")
  })
}

# Combine the list of data frames into a single large data frame
scm_all_data <- bind_rows(scm_data_list)

# -------------------------------------------------------------------------
# Part 3: Calculate Summary Statistics Across All Simulations
# -------------------------------------------------------------------------
cat("\nAggregating results and calculating summary statistics...\n")

# Define the numeric columns for which to calculate statistics
numeric_cols <- c("branch_length_D", "branch_length_I", "transitions_DtoI", 
                  "transitions_ItoD", "total_branch_length", "prop_D", 
                  "prop_I", "rate_DtoI_per_My", "rate_ItoD_per_My")

# Define the summary functions to apply
summary_funs <- list(
  mean = ~mean(.x, na.rm = TRUE),
  sd = ~sd(.x, na.rm = TRUE),
  min = ~min(.x, na.rm = TRUE),
  max = ~max(.x, na.rm = TRUE)
)

# Group by time bin and calculate statistics for all numeric columns
scm_summary <- scm_all_data %>%
  group_by(bin_id, bin_start_age, bin_end_age, midpoint_age) %>%
  summarise(
    n_simulations = n(),
    # Use across() to apply summary_funs to all columns in numeric_cols
    across(all_of(numeric_cols), summary_funs, .names = "{.col}_{.fn}"),
    .groups = "drop" # Ungroup after summarising
  )

# -------------------------------------------------------------------------
# Part 4: Data Validation and Saving
# -------------------------------------------------------------------------

# Check if proportion values are within the expected [0, 1] range
prop_cols <- c("prop_D_mean", "prop_I_mean")
for(col in prop_cols) {
  if (any(scm_summary[[col]] < 0 | scm_summary[[col]] > 1, na.rm = TRUE)) {
    out_of_range_count <- sum(scm_summary[[col]] < 0 | scm_summary[[col]] > 1, na.rm = TRUE)
    cat("Warning:", col, "has", out_of_range_count, "values outside the [0, 1] range.\n")
  }
}

# Check if mean rates are negative (which would be biologically nonsensical)
rate_cols <- c("rate_DtoI_per_My_mean", "rate_ItoD_per_My_mean")
for(col in rate_cols) {
  if (any(scm_summary[[col]] < 0, na.rm = TRUE)) {
    negative_count <- sum(scm_summary[[col]] < 0, na.rm = TRUE)
    cat("Warning:", col, "has", negative_count, "negative values.\n")
  }
}

# Save the final summary table
write_csv(scm_summary, output_file)

cat("\nSCM summary statistics have been saved to:", output_file, "\n")
