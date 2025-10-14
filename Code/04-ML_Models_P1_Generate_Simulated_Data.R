# =========================================================================
# Generate Randomized Simulated Datasets for Null Model Comparison
# Last updated: 2025-10-14
# =========================================================================

# Set the primary working directory for the project.
setwd("D:/R_work/Inflorescence")

# 1. Load Original Data
# This script uses the observed data as a template to create simulated datasets.
data_path <- "MLModel/ObData/ObData.csv"
data <- read.csv(data_path, stringsAsFactors = FALSE)

# 2. Define and Subset Variables
# Define the dependent variables (response variables).
dependent_vars <- c("H_Rank", "Indeterminate_adj")

# Define the independent environmental variables (predictors).
environmental_vars <- c(
  "Altitude", "Slope", "AMT", "AP", "MDTR", "ISO", "TS", "PS", "MTWQ",
  "PDQ", "RAD", "WIND", "AI", "ET0", "AMTd", "APd", "AMTsd", "APsd", "CSI",
  "BULK", "ECE", "OC", "PH", "SAND", "Soildiv"
)

# Combine column names for subsetting.
selected_cols <- c(dependent_vars, environmental_vars)

# Create a new dataframe containing only the selected columns.
new_data <- data[, selected_cols]

# 3. Specify Columns for Shuffling
# The environmental variables are the predictors to be randomized.
shuffle_cols <- environmental_vars

# 4. Generate and Save 100 Shuffled Datasets
# The output directory for the simulated data.
output_dir <- "MLModel/SimData/"

# Create the output directory if it does not exist.
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Get the number of rows for sampling.
n_rows <- nrow(new_data)

# Loop to generate 100 datasets with shuffled predictors.
for (i in 1:100) {
  # Set a unique seed for each iteration to ensure reproducibility.
  set.seed(i)
  
  # Create a copy of the subsetted data for this iteration.
  shuffled_data <- new_data
  
  # Shuffle each specified environmental variable column independently.
  for (col in shuffle_cols) {
    shuffled_data[[col]] <- sample(shuffled_data[[col]], size = n_rows, replace = FALSE)
  }
  
  # Generate a unique filename corresponding to the seed.
  filename <- sprintf("SimData_%03d.csv", i)
  file_path <- file.path(output_dir, filename)
  
  # Save the shuffled data to a CSV file without row names.
  write.csv(shuffled_data, file_path, row.names = FALSE)
}

cat("All simulated datasets have been successfully generated and saved to '", output_dir, "'\n")
