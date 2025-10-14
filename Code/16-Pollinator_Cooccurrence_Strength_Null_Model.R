# =========================================================================
# Calculate Standardized Effect Size (SES) of Plant-Pollinator Co-occurrence
# This SES is termed the Co-occurrence Intensity (COI).
# Last updated: 2025-10-14
# =========================================================================

# --- Load Packages ---
library(dplyr)
library(tidyr)

# Set the primary working directory for the project.
setwd("D:/R_work/Inflorescence")

# Define a default value operator (often loaded with tidyverse packages)
# This is included from your original code's dependencies.
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}


# -------------------------------------------------------------------------
# Part 1: Read All Required Data
# -------------------------------------------------------------------------
# Observed OccNum values for plants (from script 15)
observed_occnum_df <- read.csv("Potential_Pollinators/OccNum_Plants_df_Obs.csv")
# Plant distribution data
flower_dist <- read.csv("Potential_Pollinators/OCC_Plants_China.csv")
# Pollinator distribution data
pollinator_dist <- read.csv("Potential_Pollinators/OCC_Potential_Pollinators_China.csv")
# Pollinator functional group data
pollinator_types <- read.csv("Potential_Pollinators/Type_Potential_Pollinators_China.csv")


# -------------------------------------------------------------------------
# Part 2: Data Preprocessing and Preparation for Simulation
# -------------------------------------------------------------------------
# Calculate the number of grid cells occupied by each plant species (CellNum)
cellnum_data <- flower_dist %>%
  distinct(FlowerID, Grid_ID) %>%
  group_by(FlowerID) %>%
  summarise(CellNum = n(), .groups = "drop")

# --- Prepare pollinator data for rapid lookups to speed up simulations ---
# Create a list of pollinator IDs grouped by grid cell ID
pollinators_by_grid <- split(pollinator_dist$CNPollinationID, pollinator_dist$S_ID)
# Create a named vector for mapping pollinator IDs to their functional groups
pollinator_type_map <- setNames(pollinator_types$Pollination, pollinator_types$CNPollinationID)

# Pre-calculate the total number of species in each pollinator guild
total_species_counts <- pollinator_types %>%
  group_by(Pollination) %>%
  summarise(TotalSpecies = n_distinct(CNPollinationID), .groups = "drop")
# Total number of all pollinator species
total_species_all <- n_distinct(pollinator_types$CNPollinationID)

# Determine the unique CellNum values that need to be simulated
unique_cell_nums <- sort(unique(cellnum_data$CellNum))


# -------------------------------------------------------------------------
# Part 3: Run the Null Model Simulation
# -------------------------------------------------------------------------
# Initialize a dataframe to store the simulation results (mean and sd for each CellNum)
simulation_summary <- data.frame()

# Loop through each unique CellNum value
for (current_cell_num in unique_cell_nums) {
  
  # Initialize a matrix to store results for 1000 simulations for the current CellNum
  # Rows are simulation runs, columns are the different OccNum metrics
  results_matrix <- matrix(nrow = 1000, ncol = 7)
  colnames(results_matrix) <- c("All", "Beetles", "Flies", "Bees", "Wasps", "Butterflies", "Moths")
  
  # Inner loop: 1000 simulations for the current CellNum
  for (i in 1:1000) {
    # Randomly sample grid cells equal to the plant's occupancy
    selected_cells <- sample(names(pollinators_by_grid), current_cell_num, replace = FALSE)
    # Get all unique pollinator IDs found in these randomly selected cells
    simulated_pollinators <- unique(unlist(pollinators_by_grid[selected_cells]))
    
    if (length(simulated_pollinators) > 0) {
      # Get the functional types of these simulated pollinators
      simulated_types <- na.omit(pollinator_type_map[simulated_pollinators])
      
      # Calculate OccNum_All for the simulation
      results_matrix[i, "All"] <- length(simulated_pollinators) / total_species_all
      
      # Calculate OccNum for each functional guild
      type_counts <- table(simulated_types)
      results_matrix[i, "Beetles"] <- (type_counts["Beetles"] %||% 0) / (total_species_counts$TotalSpecies[total_species_counts$Pollination == "Beetles"] %||% Inf)
      results_matrix[i, "Flies"] <- (type_counts["Flies"] %||% 0) / (total_species_counts$TotalSpecies[total_species_counts$Pollination == "Flies"] %||% Inf)
      results_matrix[i, "Bees"] <- (type_counts["Bees"] %||% 0) / (total_species_counts$TotalSpecies[total_species_counts$Pollination == "Bees"] %||% Inf)
      results_matrix[i, "Wasps"] <- (type_counts["Wasps"] %||% 0) / (total_species_counts$TotalSpecies[total_species_counts$Pollination == "Wasps"] %||% Inf)
      results_matrix[i, "Butterflies"] <- (type_counts["Butterflies"] %||% 0) / (total_species_counts$TotalSpecies[total_species_counts$Pollination == "Butterflies"] %||% Inf)
      results_matrix[i, "Moths"] <- (type_counts["Moths"] %||% 0) / (total_species_counts$TotalSpecies[total_species_counts$Pollination == "Moths"] %||% Inf)
    } else {
      # If no pollinators are found in the random cells, all OccNums are 0
      results_matrix[i, ] <- 0
    }
  }
  
  # Calculate the mean and standard deviation from the 1000 simulations
  means <- colMeans(results_matrix, na.rm = TRUE)
  sds <- apply(results_matrix, 2, sd, na.rm = TRUE)
  
  # Store the results in the summary dataframe
  summary_row <- data.frame(
    CellNum = current_cell_num,
    OccNum_All_mean = means["All"], OccNum_All_sd = sds["All"],
    OccNum_Beetles_mean = means["Beetles"], OccNum_Beetles_sd = sds["Beetles"],
    OccNum_Flies_mean = means["Flies"], OccNum_Flies_sd = sds["Flies"],
    OccNum_Bees_mean = means["Bees"], OccNum_Bees_sd = sds["Bees"],
    OccNum_Wasps_mean = means["Wasps"], OccNum_Wasps_sd = sds["Wasps"],
    OccNum_Butterflies_mean = means["Butterflies"], OccNum_Butterflies_sd = sds["Butterflies"],
    OccNum_Moths_mean = means["Moths"], OccNum_Moths_sd = sds["Moths"]
  )
  simulation_summary <- rbind(simulation_summary, summary_row)
  
  cat("Simulation completed for CellNum =", current_cell_num, "\n")
}

# -------------------------------------------------------------------------
# Part 4: Merge Data and Calculate SES (COI)
# -------------------------------------------------------------------------
# Merge CellNum information into the observed data
observed_with_cellnum <- left_join(observed_occnum_df, cellnum_data, by = "FlowerID")

# Merge the simulation results (mean and sd) into the observed data
final_data <- left_join(observed_with_cellnum, simulation_summary, by = "CellNum")

# Calculate SES (Standardized Effect Size) for each OccNum metric
final_ses_df <- final_data %>%
  mutate(
    COI_All = (OccNum_All - OccNum_All_mean) / ifelse(OccNum_All_sd == 0, 1, OccNum_All_sd),
    COI_Beetles = (OccNum_Beetles - OccNum_Beetles_mean) / ifelse(OccNum_Beetles_sd == 0, 1, OccNum_Beetles_sd),
    COI_Flies = (OccNum_Flies - OccNum_Flies_mean) / ifelse(OccNum_Flies_sd == 0, 1, OccNum_Flies_sd),
    COI_Bees = (OccNum_Bees - OccNum_Bees_mean) / ifelse(OccNum_Bees_sd == 0, 1, OccNum_Bees_sd),
    COI_Wasps = (OccNum_Wasps - OccNum_Wasps_mean) / ifelse(OccNum_Wasps_sd == 0, 1, OccNum_Wasps_sd),
    COI_Butterflies = (OccNum_Butterflies - OccNum_Butterflies_mean) / ifelse(OccNum_Butterflies_sd == 0, 1, OccNum_Butterflies_sd),
    COI_Moths = (OccNum_Moths - OccNum_Moths_mean) / ifelse(OccNum_Moths_sd == 0, 1, OccNum_Moths_sd)
  ) %>%
  select(
    FlowerID, Inflorescence, 
    COI_All, COI_Beetles, COI_Flies, COI_Bees, COI_Wasps, COI_Butterflies, COI_Moths
  )

# -------------------------------------------------------------------------
# Part 5: Save the  Results
# -------------------------------------------------------------------------
output_file_path <- "Potential_Pollinators/OccNum_Plants_df_SES.csv"
write.csv(final_ses_df, file = output_file_path, row.names = FALSE)
cat(" results have been successfully saved to:\n", normalizePath(output_file_path), "\n")

# -------------------------------------------------------------------------
# Part 6: Calculate and Save Average COI per Grid Cell (MCI)
# -------------------------------------------------------------------------
# Get unique plant-grid links
plant_grid_links <- flower_dist %>%
  select(FlowerID, Grid_ID) %>%
  distinct()

# Merge the plant-level COI values to their distribution locations
merged_coi_dist_data <- left_join(plant_grid_links, 
                                  final_ses_df %>% select(FlowerID, COI_All), 
                                  by = "FlowerID")

# Calculate the Mean Co-occurrence Intensity (MCI) per grid cell
grid_average_coi <- merged_coi_dist_data %>%
  group_by(Grid_ID) %>%
  summarise(
    Num_Plant_Species = n(),
    MCI = mean(COI_All, na.rm = TRUE),
    .groups = 'drop'
  )

# Save the grid-level average COI (MCI) results
grid_coi_output_path <- "Potential_Pollinators/Grid_Average_COI_MCI.csv"
write.csv(grid_average_coi, grid_coi_output_path, row.names = FALSE)

cat("Grid-level average COI (MCI) calculation complete. Results saved to:\n", normalizePath(grid_coi_output_path), "\n")
