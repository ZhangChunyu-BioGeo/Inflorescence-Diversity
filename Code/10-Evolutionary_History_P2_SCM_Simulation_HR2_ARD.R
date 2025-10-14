# =========================================================================
# Perform 1000 SCM Simulations for the HR2_ARD Model and Summarize Results
# Last updated: 2025-10-14
# =========================================================================

# Load required packages
library(ape)
library(dplyr)
library(phytools)
library(R.utils)

# Set the primary working directory for the project.
setwd("D:/R_work/Inflorescence")

# -------------------------------------------------------------------------
# Part 1: Parameters for HR2_ARD Model and Simulation
# -------------------------------------------------------------------------

# Transition rates derived from the best-fit HR2_ARD model
RATE_D1_TO_I1 <- 0.125779788
RATE_I1_TO_D1 <- 0.069843519
RATE_D2_TO_I2 <- 0.002255834
RATE_I2_TO_D2 <- 0.0003186772
RATE_R1_TO_R2 <- 0.0278489092
RATE_R2_TO_R1 <- 0.006496064

# Seed range for simulations (defines the total number of simulations)
START_SEED <- 1
END_SEED <- 1000

# Timeout for each simulation in seconds (e.g., 6000s = 100 minutes)
# This prevents the script from hanging on a single difficult simulation.
TIMEOUT_SECONDS <- 6000

# Input files and output directories
tree_file <- "VPhyloMaker2_China_FlowerID.tre"
prob_matrix_file <- "Evolution_Model/Tip_State_Probabilities.csv"

results_dir <- "Evolution_Model"
scm_results_dir <- file.path(results_dir, "SCM_Simulations_HR2_ARD")
scm_failed_dir <- file.path(results_dir, "SCM_Simulations_HR2_ARD_Failed_Seed")
failed_seed_log_file <- file.path(results_dir, "SCM_Simulations_HR2_ARD_Failed_Seed.csv")

# Analysis parameters
NUM_BINS <- 25
TRAIT_LEVELS_OBSERVED <- c("D", "I")
TRAIT_LEVELS_HIDDEN <- c("D.R1", "I.R1", "D.R2", "I.R2")

# Create output directories
dir.create(scm_results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(scm_failed_dir, showWarnings = FALSE, recursive = TRUE) 

# -------------------------------------------------------------------------
# Part 2: Data Loading and Preprocessing
# -------------------------------------------------------------------------

phylo_tree <- read.tree(tree_file)

if (!file.exists(prob_matrix_file)) {
  stop("Error: Tip state probability matrix file not found at ", prob_matrix_file)
}
tip_probabilities_full <- read.csv(prob_matrix_file, header = TRUE, stringsAsFactors = FALSE)
colnames(tip_probabilities_full)[1] <- "Species"

# Prune tree and data to match
species_in_common <- intersect(phylo_tree$tip.label, tip_probabilities_full$Species)
pruned_tree <- keep.tip(phylo_tree, species_in_common)
tip_probabilities_pruned <- tip_probabilities_full %>% 
  filter(Species %in% species_in_common)

cat("Data preprocessing complete. Analysis is based on", length(pruned_tree$tip.label), "species.\n\n")

# Prepare the probability matrix for make.simmap (as input 'x')
tip_states_matrix <- as.matrix(tip_probabilities_pruned[, 2:5])
rownames(tip_states_matrix) <- tip_probabilities_pruned$Species
tip_states_matrix <- tip_states_matrix[pruned_tree$tip.label, ] # Ensure order matches tree
colnames(tip_states_matrix) <- TRAIT_LEVELS_HIDDEN

# Verification
if(any(abs(rowSums(tip_states_matrix) - 1) > 1e-6)) {
  warning("Row sums of the probability matrix are not exactly 1. Check for numerical precision issues.")
}

cat("Tip probability matrix prepared. Dimensions:", nrow(tip_states_matrix), "x", ncol(tip_states_matrix), "\n\n")

# -------------------------------------------------------------------------
# Part 3: Construct the Q-Matrix for the HR2_ARD Model
# -------------------------------------------------------------------------
cat("--- Constructing the HR2_ARD transition rate matrix (Q-matrix) ---\n")

q_matrix <- matrix(0, nrow = 4, ncol = 4, 
                   dimnames = list(TRAIT_LEVELS_HIDDEN, TRAIT_LEVELS_HIDDEN))

# Within-rate-category transitions (D <-> I)
q_matrix["D.R1", "I.R1"] <- RATE_D1_TO_I1
q_matrix["I.R1", "D.R1"] <- RATE_I1_TO_D1
q_matrix["D.R2", "I.R2"] <- RATE_D2_TO_I2
q_matrix["I.R2", "D.R2"] <- RATE_I2_TO_D2

# Between-rate-category transitions (R1 <-> R2)
q_matrix["D.R1", "D.R2"] <- RATE_R1_TO_R2
q_matrix["I.R1", "I.R2"] <- RATE_R1_TO_R2 # Assumes rate shift is independent of trait state
q_matrix["D.R2", "D.R1"] <- RATE_R2_TO_R1
q_matrix["I.R2", "I.R1"] <- RATE_R2_TO_R1

# Diagonals must equal the negative sum of the off-diagonal elements in each row
diag(q_matrix) <- -rowSums(q_matrix)

cat("HR2_ARD Q-matrix:\n")
print(q_matrix)

# Verification checks
if(any(is.na(q_matrix))) { stop("Q-matrix contains NA values.") }
if(any(diag(q_matrix) > 1e-10)) { stop("Diagonal elements of Q-matrix must be negative or zero.") }
if(any(abs(rowSums(q_matrix)) > 1e-10)) { stop("Row sums of Q-matrix must be zero.") }

cat("Q-matrix validation successful!\n\n")

# -------------------------------------------------------------------------
# Part 4: SCM Simulation Setup and Main Loop
# -------------------------------------------------------------------------
cat("--- Starting SCM Stochastic Simulations (HR2_ARD Model) ---\n")

# Prepare time bins for summarizing results
tree_height <- max(node.depth.edgelength(pruned_tree))
time_bins <- seq(0, tree_height, length.out = NUM_BINS + 1)
bin_midpoints <- (time_bins[-1] + time_bins[-(NUM_BINS + 1)]) / 2
node_ages <- tree_height - node.depth.edgelength(pruned_tree)

cat("Tree height:", round(tree_height, 2), "Ma\n")
cat("Time bins range from 0 to", round(tree_height, 2), "Ma\n\n")

# Create a template dataframe for failed simulations
failed_stats_template <- data.frame(
  bin_id = 1:NUM_BINS, 
  bin_start_age = time_bins[-(NUM_BINS + 1)],
  bin_end_age = time_bins[-1], 
  midpoint_age = bin_midpoints,
  branch_length_D = NA, branch_length_I = NA,
  transitions_DtoI = NA, transitions_ItoD = NA,
  total_branch_length = NA, prop_D = NA, prop_I = NA,
  rate_DtoI_per_My = NA, rate_ItoD_per_My = NA
)

# Function to handle simulation failures gracefully
handle_simulation_failure <- function(seed, reason, template_df, log_file, failed_dir) {
  log_entry <- data.frame(Timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                          Seed = seed, 
                          Reason = reason)
  write.table(log_entry, file = log_file, sep = ",", 
              append = TRUE, row.names = FALSE, 
              col.names = !file.exists(log_file))
  failed_output_csv <- file.path(failed_dir, sprintf("scm_summary_seed_%04d.csv", seed))
  write.csv(template_df, failed_output_csv, row.names = FALSE)
}

# --- Main Simulation Loop ---
success_count <- 0
failed_count <- 0
skipped_count <- 0

for (i in START_SEED:END_SEED) {
  seed_str <- sprintf("%04d", i)
  output_csv <- file.path(scm_results_dir, paste0("scm_summary_seed_", seed_str, ".csv"))
  failed_output_csv <- file.path(scm_failed_dir, paste0("scm_summary_seed_", seed_str, ".csv"))
  
  # Skip if a result (success or failure) already exists for this seed
  if (file.exists(output_csv) || file.exists(failed_output_csv)) { 
    skipped_count <- skipped_count + 1
    next
  }
  
  cat("Processing SCM simulation, Seed:", seed_str, "...\n")
  set.seed(i)
  
  simmap <- NULL
  
  tryCatch({
    # Run simulation with a timeout
    simmap_result <- withTimeout({
      make.simmap(
        tree = pruned_tree, 
        x = tip_states_matrix,  
        Q = q_matrix,
        model = "ARD", 
        nsim = 1,
        pi = "fitzjohn", # Estimate root state probabilities
        message = FALSE
      )
    }, timeout = TIMEOUT_SECONDS, onTimeout = "error")
    
    # Extract the single simmap object
    simmap <- simmap_result[[1]]
    if(is.null(simmap) || !inherits(simmap, "phylo") || is.null(simmap$maps)) {
      stop("make.simmap did not return a valid simmap object.")
    }
    
  }, TimeoutException = function(e) {
    cat("Seed", i, "timed out (>", TIMEOUT_SECONDS/60, "min). Logging failure...\n")
    handle_simulation_failure(i, "Timeout", failed_stats_template, failed_seed_log_file, scm_failed_dir)
    failed_count <<- failed_count + 1
  }, error = function(e) {
    clean_error_message <- gsub("\n", " ", e$message)
    cat("Seed", i, "failed with error:", clean_error_message, ". Logging failure...\n")
    handle_simulation_failure(i, clean_error_message, failed_stats_template, failed_seed_log_file, scm_failed_dir)
    failed_count <<- failed_count + 1
  })
  
  if(is.null(simmap)) next # Move to next iteration if simulation failed

  # -----------------------------------------------------------------------
  # Part 5: Post-processing of a single successful simulation
  # -----------------------------------------------------------------------
  bin_stats <- data.frame(
    bin_id = 1:NUM_BINS, bin_start_age = time_bins[-(NUM_BINS + 1)],
    bin_end_age = time_bins[-1], midpoint_age = bin_midpoints,
    branch_length_D = 0, branch_length_I = 0,
    transitions_DtoI = 0, transitions_ItoD = 0
  )
  
  tryCatch({
    for (j in 1:length(simmap$maps)) {
      parent_node <- simmap$edge[j, 1]
      edge_map <- simmap$maps[[j]]
      current_age <- node_ages[parent_node]
      
      if(is.null(edge_map) || length(edge_map) == 0) next
      
      for (k in 1:length(edge_map)) {
        state <- names(edge_map)[k]
        duration <- as.numeric(edge_map[k])
        
        if(is.null(state) || is.na(duration) || duration <= 0) next
        
        seg_start_age <- current_age
        seg_end_age <- current_age - duration
        
        # Determine which bins this branch segment overlaps with
        start_bin <- findInterval(seg_end_age, time_bins, rightmost.closed = TRUE)
        end_bin <- findInterval(seg_start_age, time_bins, rightmost.closed = TRUE)
        if (start_bin == 0) start_bin <- 1
        
        # Add segment duration to the correct state in each overlapping bin
        for (bin_idx in start_bin:end_bin) {
          if(bin_idx > NUM_BINS || bin_idx < 1) next
          overlap_start <- max(seg_end_age, bin_stats$bin_start_age[bin_idx])
          overlap_end <- min(seg_start_age, bin_stats$bin_end_age[bin_idx])
          overlap_duration <- overlap_end - overlap_start
          
          if (overlap_duration > 0) {
            if (state %in% c("D.R1", "D.R2")) {
              bin_stats$branch_length_D[bin_idx] <- bin_stats$branch_length_D[bin_idx] + overlap_duration
            } else if (state %in% c("I.R1", "I.R2")) {
              bin_stats$branch_length_I[bin_idx] <- bin_stats$branch_length_I[bin_idx] + overlap_duration
            }
          }
        }
        
        # Count transitions between observed states (D and I)
        if (k > 1) {
          prev_state <- names(edge_map)[k-1]
          trans_bin <- findInterval(seg_start_age, time_bins, rightmost.closed = TRUE)
          
          if (trans_bin > 0 && trans_bin <= NUM_BINS && !is.null(prev_state)) {
            observed_prev_state <- substr(prev_state, 1, 1)
            observed_state <- substr(state, 1, 1)
            
            if (observed_prev_state != observed_state) {
              if (observed_prev_state == "D" && observed_state == "I") {
                bin_stats$transitions_DtoI[trans_bin] <- bin_stats$transitions_DtoI[trans_bin] + 1
              } else if (observed_prev_state == "I" && observed_state == "D") {
                bin_stats$transitions_ItoD[trans_bin] <- bin_stats$transitions_ItoD[trans_bin] + 1
              }
            }
          }
        }
        current_age <- seg_end_age
      }
    }
    
    # Calculate final rates and proportions
    final_stats <- bin_stats %>%
      mutate(
        total_branch_length = branch_length_D + branch_length_I,
        prop_D = ifelse(total_branch_length > 1e-9, branch_length_D / total_branch_length, 0),
        prop_I = ifelse(total_branch_length > 1e-9, branch_length_I / total_branch_length, 0),
        rate_DtoI_per_My = ifelse(branch_length_D > 1e-9, transitions_DtoI / branch_length_D, 0),
        rate_ItoD_per_My = ifelse(branch_length_I > 1e-9, transitions_ItoD / branch_length_I, 0)
      )
    
    write.csv(final_stats, output_csv, row.names = FALSE)
    success_count <- success_count + 1
    
  }, error = function(e) {
    clean_error_message <- gsub("\n", " ", e$message)
    cat("Error processing results for seed", i, ":", clean_error_message, ". Logging failure...\n")
    handle_simulation_failure(i, paste("Post-processing error:", clean_error_message), failed_stats_template, failed_seed_log_file, scm_failed_dir)
    failed_count <<- failed_count + 1
  })
}

cat("\n--- Simulation Summary ---\n")
cat("Successful simulations:", success_count, "\n")
cat("Failed simulations:", failed_count, "\n")
cat("Skipped (pre-existing):", skipped_count, "\n")
cat("Total seeds processed:", START_SEED, "to", END_SEED, "\n\n")
cat("Success results saved in:", basename(scm_results_dir), "\n")
cat("Failure results saved in:", basename(scm_failed_dir), "\n")
cat("Failure log saved to:", basename(failed_seed_log_file), "\n")
