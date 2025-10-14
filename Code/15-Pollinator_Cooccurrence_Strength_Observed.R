# =========================================================================
# Calculate Observed Co-occurrence Strength between Plants and Pollinators
# Last updated: 2025-10-14
# =========================================================================

# -------------------------------------------------------------------------
# Part 1: Load Packages
# -------------------------------------------------------------------------
library(bigmemory)
library(dplyr)
library(tidyr)

# Set the primary working directory for the project.
setwd("D:/R_work/Inflorescence")

# -------------------------------------------------------------------------
# Part 2: Read Source Data
# -------------------------------------------------------------------------
# Distribution data
plant_dist <- read.csv("Potential_Pollinators/OCC_Plants_China.csv")
pollinator_dist <- read.csv("Potential_Pollinators/OCC_Potential_Pollinators_China.csv")

# Trait/Type data
plant_traits <- read.csv("Species_Trait.csv")
pollinator_types <- read.csv("Potential_Pollinators/Type_Potential_Pollinators_China.csv")

# Define the directory for big.matrix backing files
backing_dir <- "Potential_Pollinators/bigmatrix_backing/"
dir.create(backing_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------
# Part 3: Unify and Aggregate Distribution Data
# -------------------------------------------------------------------------
# Get unique species per grid cell
plant_grid_presence <- plant_dist %>%
  distinct(Grid_ID, FlowerID)

pollinator_grid_presence <- pollinator_dist %>%
  distinct(Grid_ID, CNPollinationID)

# -------------------------------------------------------------------------
# Part 4: Build and Populate the Co-occurrence Matrix
# This uses the bigmemory package to handle a potentially massive matrix
# that may not fit into RAM by storing it on disk.
# -------------------------------------------------------------------------
unique_plants <- unique(plant_grid_presence$FlowerID)
unique_pollinators <- unique(pollinator_grid_presence$CNPollinationID)

cat("Initializing a", length(unique_plants), "x", length(unique_pollinators), "file-backed co-occurrence matrix...\n")

# Create the file-backed big.matrix
co_occurrence_matrix <- bigmemory::big.matrix(
  nrow = length(unique_plants),
  ncol = length(unique_pollinators),
  type = "integer",
  init = 0,
  backingpath = backing_dir,
  backingfile = "co_occurrence_matrix.bin",
  descriptorfile = "co_occurrence_matrix.desc"
)

# Create a mapping from species ID to matrix index for fast lookups
plant_id_map <- setNames(seq_along(unique_plants), unique_plants)
pollinator_id_map <- setNames(seq_along(unique_pollinators), unique_pollinators)

# Group plant data by Grid_ID for efficient iteration
plant_grid_groups <- split(plant_grid_presence$FlowerID, plant_grid_presence$Grid_ID)
all_grid_ids <- union(names(plant_grid_groups), unique(pollinator_grid_presence$Grid_ID))
total_grids <- length(all_grid_ids)

cat("Processing", total_grids, "grid cells to count co-occurrences...\n")
progress_counter <- 0

# Iterate through each grid cell to populate the matrix
for (grid_id in all_grid_ids) {
  plants_in_grid <- plant_grid_groups[[as.character(grid_id)]]
  pollinators_in_grid <- pollinator_grid_presence$CNPollinationID[pollinator_grid_presence$Grid_ID == grid_id]
  
  # Only process if both plants and pollinators are present
  if (length(plants_in_grid) > 0 && length(pollinators_in_grid) > 0) {
    plant_indices <- plant_id_map[plants_in_grid]
    pollinator_indices <- pollinator_id_map[pollinators_in_grid]
    
    # Increment the counter for each co-occurring pair
    for (f_idx in plant_indices) {
      for (p_idx in pollinator_indices) {
        co_occurrence_matrix[f_idx, p_idx] <- co_occurrence_matrix[f_idx, p_idx] + 1
      }
    }
  }
  
  # Print progress
  progress_counter <- progress_counter + 1
  if (progress_counter %% 1000 == 0 || progress_counter == total_grids) {
    cat("Processed", progress_counter, "/", total_grids, "grids (", 
        round(progress_counter / total_grids * 100, 1), "%)\n")
  }
}

# -------------------------------------------------------------------------
# Part 5: Create Species-Pair Dataframe and Calculate OccNum
# -------------------------------------------------------------------------
cat("\nConstructing species-pair dataframe to calculate final OccNum...\n")

# Create a dataframe of all possible plant-pollinator pairs
all_pairs_df <- expand.grid(FlowerID = unique_plants,
                            pollinator_id = unique_pollinators,
                            stringsAsFactors = FALSE)

# Extract co-occurrence info and convert to binary (1 if co-occurred at least once, 0 otherwise)
all_pairs_df$co_occurred <- mapply(function(f_id, p_id) {
  f_idx <- plant_id_map[[f_id]]
  p_idx <- pollinator_id_map[[p_id]]
  return(ifelse(co_occurrence_matrix[f_idx, p_idx] > 0, 1, 0))
}, all_pairs_df$FlowerID, all_pairs_df$pollinator_id)

# Merge plant and pollinator trait/type data
all_pairs_df <- all_pairs_df %>%
  left_join(plant_traits, by = "FlowerID") %>%
  left_join(pollinator_types, by = c("pollinator_id" = "CNPollinationID")) %>%
  # Filter out pairs where trait or type info is missing
  filter(!is.na(Trait) & !is.na(Pollination))

# Calculate the final plant-level OccNum (co-occurrence strength)
# This is the proportion of potential pollinators that a plant species co-occurs with.
cat("Calculating plant-level OccNum for all pollinator guilds...\n")
OccNum_Plants_df_Obs <- all_pairs_df %>%
  group_by(FlowerID) %>%
  summarise(
    Trait = first(Trait),
    OccNum_All = mean(co_occurred, na.rm = TRUE),
    OccNum_Beetles = mean(co_occurred[Pollination == "Beetles"], na.rm = TRUE),
    OccNum_Flies = mean(co_occurred[Pollination == "Flies"], na.rm = TRUE),
    OccNum_Bees = mean(co_occurred[Pollination == "Bees"], na.rm = TRUE),
    OccNum_Wasps = mean(co_occurred[Pollination == "Wasps"], na.rm = TRUE),
    OccNum_Butterflies = mean(co_occurred[Pollination == "Butterflies"], na.rm = TRUE),
    OccNum_Moths = mean(co_occurred[Pollination == "Moths"], na.rm = TRUE),
    .groups = 'drop'
  )

# Replace any NaN (from dividing by zero, e.g., no beetles in data) with 0
OccNum_Plants_df_Obs <- OccNum_Plants_df_Obs %>%
  mutate(across(starts_with("OccNum_"), ~ifelse(is.nan(.), 0, .)))

# -------------------------------------------------------------------------
# Part 6: Save the  Results
# -------------------------------------------------------------------------
output_file_path <- "Potential_Pollinators/OccNum_Plants_df_Obs.csv"
write.csv(OccNum_Plants_df_Obs, file = output_file_path, row.names = FALSE)
cat("\nObserved co-occurrence strength results successfully saved to:\n", 
    normalizePath(output_file_path), "\n")
