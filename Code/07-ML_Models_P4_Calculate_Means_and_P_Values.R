# =========================================================================
# Calculate P-Values and Means from Simulated and Observed Model Results
# Last updated: 2025-10-14
# =========================================================================

# Set the primary working directory for the project.
setwd("D:/R_work/Inflorescence")

# -------------------------------------------------------------------------
# Part 1: Path and Basic Configuration
# -------------------------------------------------------------------------
obmodel_path <- "MLModel/ObModel/"
simmodel_path <- "MLModel/SimModel/"
result_path <- "MLModel/Result/"

# List of files that require column-wise min-max normalization before t-test
# Only the GBM model's importance scores need this special handling.
normalize_csvs <- c("GBM_IMP.csv")

# Create the Result directory if it does not exist
if (!dir.exists(result_path)) {
    dir.create(result_path, recursive = TRUE)
}

# Generate subfolder names (from "001" to "100")
subfolders <- sprintf("%03d", 1:100)

# Get the list of CSV files from the first ObModel subfolder to serve as a template
first_subfolder <- file.path(obmodel_path, subfolders[1])
csv_files <- list.files(first_subfolder, pattern = "\\.csv$", full.names = FALSE)

# -------------------------------------------------------------------------
# Part 2: Min-Max Normalization Helper Function
# -------------------------------------------------------------------------
# Function to perform column-wise min-max normalization on a single dataframe
min_max_normalize_df <- function(df) {
    if (is.null(df)) return(NULL)
    for (j in seq_along(df)) {
        col_values <- df[[j]]
        min_val <- min(col_values, na.rm = TRUE)
        max_val <- max(col_values, na.rm = TRUE)
        # If max equals min, the column has constant values or is all NA; keep it as is.
        if (!is.na(min_val) && !is.na(max_val) && (max_val != min_val)) {
            df[[j]] <- (col_values - min_val) / (max_val - min_val)
        }
    }
    return(df)
}

# -------------------------------------------------------------------------
# Part 3: Main Loop to Process Each CSV File Type
# -------------------------------------------------------------------------
for (csv_file in csv_files) {
    cat("Processing CSV file:", csv_file, "\n")
    
    # Initialize lists to store data from ObModel and SimModel runs
    ob_data_list <- list()
    sim_data_list <- list()
    
    # --- 3.1 Load data from all subfolders ---
    for (subfolder in subfolders) {
        ob_csv_path <- file.path(obmodel_path, subfolder, csv_file)
        sim_csv_path <- file.path(simmodel_path, subfolder, csv_file)
        
        # Read ObModel CSV, using the first column as row names
        if (file.exists(ob_csv_path)) {
            ob_data <- try(
                read.csv(ob_csv_path, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE),
                silent = TRUE
            )
            if (inherits(ob_data, "try-error")) {
                warning(paste("Error reading file:", ob_csv_path))
                ob_data <- NULL
            }
            ob_data_list[[subfolder]] <- ob_data
        } else {
            warning(paste("File not found:", ob_csv_path))
            ob_data_list[[subfolder]] <- NULL
        }
        
        # Read SimModel CSV, using the first column as row names
        if (file.exists(sim_csv_path)) {
            sim_data <- try(
                read.csv(sim_csv_path, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE),
                silent = TRUE
            )
            if (inherits(sim_data, "try-error")) {
                warning(paste("Error reading file:", sim_csv_path))
                sim_data <- NULL
            }
            sim_data_list[[subfolder]] <- sim_data
        } else {
            warning(paste("File not found:", sim_csv_path))
            sim_data_list[[subfolder]] <- NULL
        }
    }
    
    # --- 3.2 Calculate and save the mean of ObModel results ---
    valid_ob_data <- ob_data_list[!sapply(ob_data_list, is.null)]
    if (length(valid_ob_data) == 0) {
        warning(paste("No valid ObModel data found for", csv_file))
        next
    }
    
    # Ensure row and column names are consistent across all files
    reference_ob <- valid_ob_data[[1]]
    for (i in seq_along(valid_ob_data)) {
        if (!identical(rownames(reference_ob), rownames(valid_ob_data[[i]])) ||
            !identical(colnames(reference_ob), colnames(valid_ob_data[[i]]))) {
            stop(paste("Row/column names are inconsistent for file", csv_file, "in subfolders."))
        }
    }
    
    # Convert list of dataframes to a 3D array (rows, cols, subfolders)
    ob_matrices <- lapply(valid_ob_data, as.matrix)
    ob_array <- simplify2array(ob_matrices)
    
    # Calculate the mean across the 3rd dimension (subfolders), ignoring NAs
    ob_mean <- apply(ob_array, c(1, 2), function(x) {
        if (all(is.na(x))) NA else mean(x, na.rm = TRUE)
    })
    
    # Convert back to a dataframe with original row/column names
    ob_mean_df <- as.data.frame(ob_mean)
    
    # Save the mean result CSV
    mean_csv_path <- file.path(result_path, csv_file)
    write.csv(ob_mean_df, mean_csv_path, row.names = TRUE, na = "")
    
    # --- 3.3 Perform t-tests between ObModel and SimModel results ---
    valid_sim_data <- sim_data_list[!sapply(sim_data_list, is.null)]
    if (length(valid_sim_data) == 0) {
        warning(paste("No valid SimModel data found for", csv_file))
        next
    }
    
    # Ensure consistency for SimModel data as well
    reference_sim <- valid_sim_data[[1]]
    
    # If the current file requires normalization, apply it to both ObModel and SimModel data lists
    if (csv_file %in% normalize_csvs) {
        norm_ob_data_list <- lapply(ob_data_list, min_max_normalize_df)
        norm_sim_data_list <- lapply(sim_data_list, min_max_normalize_df)
    } else {
        # Otherwise, use the original data for the t-test
        norm_ob_data_list <- ob_data_list
        norm_sim_data_list <- sim_data_list
    }
    
    # Convert the (potentially normalized) lists to 3D arrays for t-testing
    norm_valid_ob_data <- norm_ob_data_list[!sapply(norm_ob_data_list, is.null)]
    ob_array_for_test <- simplify2array(lapply(norm_valid_ob_data, as.matrix))
    
    norm_valid_sim_data <- norm_sim_data_list[!sapply(norm_sim_data_list, is.null)]
    sim_array_for_test <- simplify2array(lapply(norm_valid_sim_data, as.matrix))
    
    # Initialize a matrix to store p-values
    p_values <- matrix(NA, nrow = nrow(reference_ob), ncol = ncol(reference_ob),
                       dimnames = list(rownames(reference_ob), colnames(reference_ob)))
    
    # --- 3.4 Perform t-test for each corresponding cell ---
    for (r in 1:nrow(reference_ob)) {
        for (c in 1:ncol(reference_ob)) {
            ob_values <- ob_array_for_test[r, c, ]
            sim_values <- sim_array_for_test[r, c, ]
            
            # Remove NAs
            ob_values_clean <- ob_values[!is.na(ob_values)]
            sim_values_clean <- sim_values[!is.na(sim_values)]
            
            # Perform t-test only if there are at least 2 data points in each group
            if (length(ob_values_clean) >= 2 && length(sim_values_clean) >= 2) {
                test_result <- try(t.test(ob_values_clean, sim_values_clean), silent = TRUE)
                if (!inherits(test_result, "try-error")) {
                    p_values[r, c] <- test_result$p.value
                }
            }
        }
    }
    
    # Convert the p-value matrix to a dataframe
    p_values_df <- as.data.frame(p_values)
    
    # Generate the output filename with a "_p" suffix
    p_csv_name <- sub("\\.csv$", "_p.csv", csv_file)
    p_csv_path <- file.path(result_path, p_csv_name)
    write.csv(p_values_df, p_csv_path, row.names = TRUE, na = "")
}

cat("All files have been processed.\n")

