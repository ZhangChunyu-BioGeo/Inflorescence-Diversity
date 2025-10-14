# =========================================================================
# Fit and Compare Multiple Evolutionary Models (10 models including 
# standard, hidden-state, and variable-rate)
# Last updated: 2025-10-14
# =========================================================================

# Load R packages
library(ape)
library(corHMM)
library(dplyr)

# Set the primary working directory for the project.
setwd("D:/R_work/Inflorescence")

# -------------------------------------------------------------------------
# Part 1: Global Parameters
# -------------------------------------------------------------------------
tree_file <- "VPhyloMaker2_China_FlowerID.tre"
trait_file <- "Species_Trait.csv"

# Directory for storing all evolution-related results
results_dir <- "Evolution_Model"
# Subdirectory specifically for the fitted model objects
model_results_dir <- file.path(results_dir, "Model_Fits")

# Define the states/levels of the categorical trait
TRAIT_LEVELS <- c("D", "I")

# Optional: Manually specify model names here to skip them during the run
# Example: MODELS_TO_SKIP <- c("ER_H2", "Var_R1_ARD_H2")
MODELS_TO_SKIP <- c()

# Create output directories if they don't exist
dir.create(model_results_dir, showWarnings = FALSE, recursive = TRUE)

# The script will fit and select from the following ten models:
# - ER, ARD, Dir_DtoI, Dir_ItoD,
# - ER_H2, ARD_H2, Dir_DtoI_H2, Dir_ItoD_H2,
# - Var_R1_ARD_H2, Var_R1_ER_H2
if (length(MODELS_TO_SKIP) > 0) {
    cat("\nBased on settings, the following models will be skipped:\n")
    cat("  - ", paste(MODELS_TO_SKIP, collapse = ", "), "\n")
}

# -------------------------------------------------------------------------
# Part 2: Data Loading and Preprocessing
# -------------------------------------------------------------------------
phylo_tree <- read.tree(tree_file)
trait_data <- read.csv(trait_file, header = TRUE, stringsAsFactors = FALSE)

# Standardize column names
colnames(trait_data)[colnames(trait_data) == "FlowerID"] <- "Species"
trait_data$Trait <- factor(trait_data$Trait, levels = TRAIT_LEVELS)

# Prune tree and data to include only common species
species_in_common <- intersect(phylo_tree$tip.label, trait_data$Species)
pruned_tree <- keep.tip(phylo_tree, species_in_common)
trait_data_pruned <- trait_data %>% filter(Species %in% species_in_common)

# Format data for corHMM (species name and numeric trait)
data_corHMM <- trait_data_pruned[, c("Species", "Trait")]
data_corHMM$Trait <- as.numeric(data_corHMM$Trait)

cat("Data preprocessing complete. Analysis will be based on", length(pruned_tree$tip.label), "species.\n\n")

# -------------------------------------------------------------------------
# Part 3: Model Definition and Fitting
# -------------------------------------------------------------------------

# 1. Define four base rate matrices (for models with 1 rate category)
base_matrices <- list(
  ER = corHMM::getStateMat4Dat(data_corHMM, "ER")$rate.mat,
  ARD = corHMM::getStateMat4Dat(data_corHMM, "ARD")$rate.mat,
  Dir_DtoI = { mat <- corHMM::getStateMat4Dat(data_corHMM, "ARD")$rate.mat; mat[2, 1] <- 0; mat },
  Dir_ItoD = { mat <- corHMM::getStateMat4Dat(data_corHMM, "ARD")$rate.mat; mat[1, 2] <- 0; mat }
)

# 2. Generate a list of all final Q matrices for all 10 models
full_rate_matrices <- list()
for (base_model in names(base_matrices)) {
  # Add the standard (1-rate) model
  full_rate_matrices[[base_model]] <- base_matrices[[base_model]]
  
  # Add the hidden-state (2-rate) model
  model_name_h2 <- paste0(base_model, "_H2")
  RateCatMat <- corHMM::getRateCatMat(2) # Matrix for transitions between rate categories
  full_rate_matrices[[model_name_h2]] <- corHMM::getFullMat(list(base_matrices[[base_model]], base_matrices[[base_model]]), RateCatMat)
}

# 2b. Manually construct two additional "variable-rate" hidden-state models
mat_var_r1_ard <- matrix(c(0,2,5,0, 1,0,0,5, 4,0,0,3, 0,4,3,0), 4, 4)
full_rate_matrices[["Var_R1_ARD_H2"]] <- mat_var_r1_ard

mat_var_r1_er <- matrix(c(0,1,5,0, 1,0,0,5, 4,0,0,3, 0,4,2,0), 4, 4)
full_rate_matrices[["Var_R1_ER_H2"]] <- mat_var_r1_er

cat("Successfully defined rate matrices for all 10 models.\n\n")

# 3. Loop to run and save all models
model_fits <- list()
new_models_fitted <- 0
existing_models_loaded <- 0
manual_skips <- 0

for (model_name in names(full_rate_matrices)) {
  
  # Check if this model should be skipped manually
  if (model_name %in% MODELS_TO_SKIP) {
    cat("Manually skipping model:", model_name, "...\n")
    manual_skips <- manual_skips + 1
    next # Proceed to the next iteration
  }
  
  rate_mat_to_use <- full_rate_matrices[[model_name]]
  num_rate_cats <- ifelse(grepl("_H2", model_name), 2, 1)
  
  # Define paths for saving/loading results
  model_dir <- file.path(model_results_dir, paste0("model_", model_name))
  dir.create(model_dir, showWarnings = FALSE)
  rds_file <- file.path(model_dir, paste0(model_name, "_fit.rds"))
  
  if (file.exists(rds_file)) {
    cat("Found existing result for model '", model_name, "', loading...\n", sep="")
    model_fits[[model_name]] <- readRDS(rds_file)
    existing_models_loaded <- existing_models_loaded + 1
  } else {
    cat("Running ", model_name, " model (rate.cat = ", num_rate_cats, ")...\n", sep="")
    start_time <- Sys.time()
    
    fit <- corHMM(
      phy = pruned_tree, 
      data = data_corHMM, 
      rate.cat = num_rate_cats, 
      rate.mat = rate_mat_to_use,
      root.p = "marginal",
      nstarts = 10,       # Number of random starting points to find the best likelihood
      n.cores = 1,        # Set to >1 for parallel processing if supported
      upper.bound = 10    # Upper bound for rate estimates
    )
    
    end_time <- Sys.time()
    runtime <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)
    
    # Save the fitted model object and a text summary
    saveRDS(fit, rds_file)
    sink(file.path(model_dir, paste0(model_name, "_summary.txt"))); print(fit); sink()
    
    model_fits[[model_name]] <- fit
    new_models_fitted <- new_models_fitted + 1
    cat(model_name, " model finished and saved. Runtime:", runtime, "minutes\n")
  }
}

# -------------------------------------------------------------------------
# Part 4: Model Comparison and Summary
# -------------------------------------------------------------------------

# Create a model comparison table based on AICc
model_comparison <- data.frame(
  Model = names(model_fits),
  LogLik = sapply(model_fits, '[[', 'loglik'),
  AICc = sapply(model_fits, '[[', 'AICc'),
  N_params = sapply(model_fits, function(x) length(x$solution[x$solution > 0])),
  stringsAsFactors = FALSE
) %>% arrange(AICc)

# Calculate delta AICc (difference from the best model)
model_comparison$Delta_AICc <- model_comparison$AICc - min(model_comparison$AICc)

# Save the summary table to a CSV file
summary_file_path <- file.path(results_dir, "Model_Selection_Summary.csv")
write.csv(model_comparison, summary_file_path, row.names = FALSE)

cat("\nModel fitting and comparison complete.\n")
cat("Summary table saved to:", summary_file_path, "\n")
cat(" - New models fitted:", new_models_fitted, "\n")
cat(" - Existing models loaded:", existing_models_loaded, "\n")
cat(" - Models manually skipped:", manual_skips, "\n")
