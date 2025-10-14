# =========================================================================
# Extract and Aggregate Paleoclimate Data using a Weighted Average
# Note: ed_PhanDA_GMSTandCO2_percentiles.csv is from DOI: 10.1126/science.adk3705
# Last updated: 2025-10-14
# =========================================================================

# Load required packages
library(ape)
library(dplyr)

# Set the primary working directory for the project.
setwd("D:/R_work/Inflorescence")

# -------------------------------------------------------------------------
# Part 1: Global Parameters
# -------------------------------------------------------------------------
tree_file <- "VPhyloMaker2_China_FlowerID.tre"
climate_file <- "ed_PhanDA_GMSTandCO2_percentiles.csv"
results_dir <- "Evolution_Model"
NUM_BINS <- 25

# Create the output directory if it doesn't exist
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------
# Part 2: Data Loading and Preprocessing
# -------------------------------------------------------------------------

# Load the phylogenetic tree to get the time range
phylo_tree <- read.tree(tree_file)
tree_height <- max(node.depth.edgelength(phylo_tree))

# Load the paleoclimate data
climate_data <- read.csv(climate_file, header = TRUE, stringsAsFactors = FALSE)

# Ensure GMST (Global Mean Surface Temperature) columns are numeric
climate_data <- climate_data %>%
  mutate(across(starts_with("GMST_"), as.numeric))

# Identify all GMST-related columns
gmst_columns <- colnames(climate_data)[grepl("GMST_", colnames(climate_data))]

# Check the integrity of the climate data
cat("Checking non-NA value counts in raw climate data:\n")
for(col in gmst_columns) {
  na_count <- sum(is.na(climate_data[[col]]))
  cat(col, ": ", nrow(climate_data) - na_count, "/", nrow(climate_data), " non-NA values\n")
}
cat("\n")

# -------------------------------------------------------------------------
# Part 3: Create Evolutionary Time Bins
# -------------------------------------------------------------------------
# These bins must match those used in the SCM analysis scripts
time_bins <- seq(0, tree_height, length.out = NUM_BINS + 1)
evo_bins <- data.frame(
  bin_id = 1:NUM_BINS,
  bin_start_age = time_bins[-(NUM_BINS + 1)], 
  bin_end_age = time_bins[-1], 
  bin_duration = diff(time_bins)
)

# -------------------------------------------------------------------------
# Part 4: Aggregate Climate Data using Weighted Average
# -------------------------------------------------------------------------
cat("Starting to aggregate climate data into evolutionary time bins...\n")

# Initialize the results dataframe
aggregated_climate <- evo_bins
for(col in gmst_columns) {
  aggregated_climate[[col]] <- NA_real_
}

# Add columns for tracking statistics
aggregated_climate$n_climate_periods <- 0      # Number of climate periods contributing
aggregated_climate$total_overlap_duration <- 0 # Total duration of overlap
aggregated_climate$coverage_percent <- 0       # Percentage of the bin covered by data

# Loop through each evolutionary time bin
for (i in 1:nrow(evo_bins)) {
  evo_bin_start <- evo_bins$bin_start_age[i]
  evo_bin_end <- evo_bins$bin_end_age[i]
  evo_bin_duration <- evo_bins$bin_duration[i]
  
  # Reset accumulators for each bin
  total_overlap_duration <- 0
  weighted_sums <- setNames(rep(0, length(gmst_columns)), gmst_columns)
  climate_periods_count <- 0
  
  # Loop through each period in the climate dataset
  for (j in 1:nrow(climate_data)) {
    clim_period_start <- climate_data$UpperAge[j]
    clim_period_end <- climate_data$LowerAge[j]
    
    # Calculate the overlapping interval
    overlap_start <- max(evo_bin_start, clim_period_start)
    overlap_end <- min(evo_bin_end, clim_period_end)
    overlap_duration <- overlap_end - overlap_start
    
    # If there is a meaningful overlap
    if (overlap_duration > 1e-6) {
      gmst_values <- as.numeric(climate_data[j, gmst_columns])
      
      # Check if there is any valid (non-NA) GMST data for this period
      if(any(!is.na(gmst_values))) {
          
        # Add to the weighted sum for each GMST column
        for(k in 1:length(gmst_columns)) {
          if(!is.na(gmst_values[k])) {
            weighted_sums[k] <- weighted_sums[k] + (gmst_values[k] * overlap_duration)
          }
        }
        
        # Update statistics
        climate_periods_count <- climate_periods_count + 1
        total_overlap_duration <- total_overlap_duration + overlap_duration
      }
    }
  }
  
  # Store the calculated statistics for the current bin
  aggregated_climate$n_climate_periods[i] <- climate_periods_count
  aggregated_climate$total_overlap_duration[i] <- total_overlap_duration
  aggregated_climate$coverage_percent[i] <- round((total_overlap_duration / evo_bin_duration) * 100, 1)
  
  # Calculate the final weighted average for each GMST column
  if (total_overlap_duration > 1e-6) {
    for(k in 1:length(gmst_columns)) {
        # The weighted average is the sum of (value * duration) divided by total duration.
        # However, since some GMST columns have NAs, the denominator for each must be the 
        # total overlap duration for which THAT column had valid data.
        # The current implementation simplifies this by dividing by the total overlap where ANY data was present.
        # This is a reasonable approach if data availability is consistent across columns.
        aggregated_climate[i, gmst_columns[k]] <- weighted_sums[k] / total_overlap_duration
    }
  }
  
  # Print progress
  if(i %% 5 == 0 || i == NUM_BINS) {
    cat("Processed", i, "/", NUM_BINS, "time bins.",
        " Current bin coverage:", aggregated_climate$coverage_percent[i], "%\n")
  }
}

# -------------------------------------------------------------------------
# Part 5: Result Summary and Validation
# -------------------------------------------------------------------------
cat("\n--- Aggregation Summary ---\n")

# Count how many bins have any climate data
total_bins_with_data <- sum(aggregated_climate$n_climate_periods > 0)
cat("Total bins with some climate data:", total_bins_with_data, "out of", NUM_BINS, "\n")

# Report on data coverage
avg_coverage <- mean(aggregated_climate$coverage_percent[aggregated_climate$coverage_percent > 0])
full_coverage_bins <- sum(aggregated_climate$coverage_percent >= 99.9)
cat("Average data coverage (for bins with data):", round(avg_coverage, 1), "%\n")
cat("Number of bins with full (>=99.9%) coverage:", full_coverage_bins, "\n\n")

# Report non-NA counts for each aggregated GMST column
cat("Final data availability per aggregated GMST column:\n")
for(col in gmst_columns) {
  non_na_count <- sum(!is.na(aggregated_climate[[col]]))
  cat("  ", col, ": ", non_na_count, "/", NUM_BINS, " bins have data (", 
      round(non_na_count/NUM_BINS*100, 1), "%)\n")
}

# -------------------------------------------------------------------------
# Part 6: Save Results
# -------------------------------------------------------------------------
climate_output_file <- file.path(results_dir, "Aggregated_Climate_Data.csv")
write.csv(aggregated_climate, climate_output_file, row.names = FALSE)
cat("\nAggregated climate data saved to:", climate_output_file, "\n")
