# =========================================================================
# Process Simulated Datasets in a Loop and Save Model Results
# Last updated: 2025-10-14
# =========================================================================

# Load necessary packages
library(caret)
library(randomForest)
library(gbm)
library(e1071)
library(nnet)
library(iml)
library(DALEX)
library(kernlab)
library(dplyr)
library(MASS)
library(purrr)
library(glmnet)
library(fastshap)

# Set the primary working directory for the project.
setwd("D:/R_work/Inflorescence")

# Define the range of seeds for the loop
seed_start <- 1
seed_end <- 100

# Define input data directory and the root path for output results
input_data_dir <- "MLModel/SimData/"
output_root_path <- "MLModel/SimModel/"

# Define the filename template for simulated data
sim_data_template <- "SimData_%03d.csv"

# Define columns to be normalized
columns_to_normalize <- c(
  "H_Rank", "Indeterminate_adj",
  "Altitude", "Slope", "AMT", "AP", "MDTR", "ISO", "TS", "PS", "MTWQ",
  "PDQ", "RAD", "WIND", "AI", "ET0", "AMTd", "APd", "AMTsd", "APsd", "CSI",
  "BULK", "ECE", "OC", "PH", "SAND", "Soildiv"
)

# Define normalization function (Min-Max scaling)
normalize <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

# Define model list and related settings
models <- c('RF', 'GBM', 'SVR', 'KNN', 'GLM', 'NN')

model_list <- list(
  RF = "rf",
  GBM = "gbm",
  SVR = "svmRadial",
  KNN = "knn",
  GLM = "glm",
  NN = "nnet"
)

train_control <- trainControl(
  method = "cv",
  number = 5,
  search = "random",
  allowParallel = FALSE
)

# Define environmental predictor variables
env_vars <- c(
  "Altitude", "Slope", "AMT", "AP", "MDTR", "ISO", "TS", "PS", "MTWQ",
  "PDQ", "RAD", "WIND", "AI", "ET0", "AMTd", "APd", "AMTsd", "APsd", "CSI",
  "BULK", "ECE", "OC", "PH", "SAND", "Soildiv"
)

# Define target (response) variable column names
target_vars <- c("H_Rank", "Indeterminate_adj")

# Define a function to save the results
save_results <- function(feature_importance_IMP, feature_importance_SHAP, Model_R2, Model_RMSE, Model_MSE, save_path, models) {
  
  # Internal function to save feature importance data
  save_feature_importance <- function(model_name, metric, data_list, save_path) {
    imp_data <- data_list[[model_name]]
    if (all(is.na(imp_data))) {
      cat("    All ", metric, " data for model ", model_name, " are NA. Skipping save.\n")
      return(NULL)
    }
    
    imp_df <- data.frame(
      Feature = rownames(imp_data),
      imp_data,
      stringsAsFactors = FALSE
    )
    
    filename <- paste0(save_path, model_name, "_", metric, ".csv")
    write.csv(imp_df, file = filename, row.names = FALSE)
    cat("    ", metric, " data for model ", model_name, " saved to ", filename, "\n")
  }
  
  # Iterate over models to save IMP and SHAP data
  for (model in models) {
    save_feature_importance(model, "IMP", feature_importance_IMP, save_path)
    save_feature_importance(model, "SHAP", feature_importance_SHAP, save_path)
  }
  
  # Save performance matrices
  write.csv(Model_R2, file = paste0(save_path, "Model_R2.csv"), row.names = TRUE)
  write.csv(Model_RMSE, file = paste0(save_path, "Model_RMSE.csv"), row.names = TRUE)
  write.csv(Model_MSE, file = paste0(save_path, "Model_MSE.csv"), row.names = TRUE)
}

# ---------------------------- MAIN LOOP ----------------------------
# Loop through each seed and its corresponding data file
for (seed in seed_start:seed_end) {
  
  seed_str <- sprintf("%03d", seed)
  cat("=== Running Seed:", seed_str, "===\n")
  
  # Construct the path to the current data file
  current_data_file <- file.path(input_data_dir, sprintf(sim_data_template, seed))
  
  if (!file.exists(current_data_file)) {
    cat("  Warning: File ", current_data_file, " does not exist. Skipping this seed.\n\n")
    next
  }
  
  # Set the seed for reproducibility for this iteration
  MySeed <- as.integer(seed_str)
  set.seed(MySeed)
  
  # Read the current data file
  df <- read.csv(file = current_data_file, header = TRUE)
  
  # Normalize the specified columns
  df[columns_to_normalize] <- lapply(df[columns_to_normalize], normalize)
  
  # Initialize lists and matrices to store results for this seed
  feature_importance_IMP <- list()
  feature_importance_SHAP <- list()
  
  init_matrix <- function() {
    matrix(NA, nrow = length(env_vars), ncol = length(target_vars),
           dimnames = list(env_vars, target_vars))
  }
  
  for (model in models) {
    feature_importance_IMP[[model]] <- init_matrix()
    feature_importance_SHAP[[model]] <- init_matrix()
  }
  
  init_performance_matrix <- function() {
    matrix(NA, nrow = length(models), ncol = length(target_vars),
           dimnames = list(models, target_vars))
  }
  
  Model_R2 <- init_performance_matrix()
  Model_RMSE <- init_performance_matrix()
  Model_MSE <- init_performance_matrix()
  
  # --- Function Definitions for Metrics ---
  compute_shap <- function(model_fit, train_data, env_vars, target, nsim = 100) {
    pred_fun <- function(model, newdata) {
      predict(model, newdata = newdata)
    }
    
    shap_values <- tryCatch({
      fastshap::explain(
        object = model_fit$finalModel,
        X = train_data[, env_vars],
        pred_wrapper = pred_fun,
        nsim = nsim
      )
    }, error = function(e) {
      cat("    SHAP calculation error:", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (is.null(shap_values)) return(NULL)
    
    shap_imp <- colMeans(abs(shap_values))
    shap_imp <- shap_imp[env_vars]
    
    return(shap_imp)
  }
  
  compute_feature_imp_iml <- function(model_fit, train_data, env_vars, target) {
    predictor <- tryCatch({
      Predictor$new(model_fit, data = train_data[, env_vars], y = train_data[[target]])
    }, error = function(e) {
      cat("    Predictor object creation error:", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (is.null(predictor)) return(NULL)
    
    feature_imp <- tryCatch({
      FeatureImp$new(predictor, loss = "rmse")
    }, error = function(e) {
      cat("    FeatureImp calculation error:", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (is.null(feature_imp)) return(NULL)
    
    imp_df <- feature_imp$results
    imp_df <- imp_df[imp_df$feature %in% env_vars, ]
    
    imp_vec <- setNames(imp_df$importance, imp_df$feature)
    imp_vec_complete <- setNames(rep(NA, length(env_vars)), env_vars)
    imp_vec_complete[names(imp_vec)] <- imp_vec
    
    return(imp_vec_complete)
  }

  # --- Process Each Target Variable ---
  for (target in target_vars) {
    cat("  Processing target variable:", target, "\n")
    flush.console()
    
    data <- df[, c(env_vars, target)]
    data <- na.omit(data)
    
    set.seed(MySeed)
    train_index <- createDataPartition(data[[target]], p = 0.7, list = FALSE)
    train_data <- data[train_index, ]
    test_data <- data[-train_index, ]
    
    train_data <- na.omit(train_data)
    test_data <- na.omit(test_data)
    
    # --- Iterate Through Each Model ---
    for (model in models) {
      method <- model_list[[model]]
      cat("    Training", model, "model...\n")
      flush.console()
      
      model_seed <- MySeed + as.numeric(as.factor(target)) * 100 + as.numeric(as.factor(model))
      set.seed(model_seed)
      
      model_fit <- tryCatch({
        train_formula <- as.formula(paste(target, "~", paste(env_vars, collapse = "+")))
        if (model == "NN") {
          caret::train(train_formula, data = train_data, method = method, metric = "RMSE",
                       trControl = train_control, tuneLength = 10, linout = TRUE, trace = FALSE)
        } else if (model == "GLM") {
          caret::train(train_formula, data = train_data, method = method, metric = "RMSE",
                       trControl = train_control, tuneLength = 10, family = gaussian())
        } else {
          caret::train(train_formula, data = train_data, method = method, metric = "RMSE",
                       trControl = train_control, tuneLength = 10, verbose = FALSE)
        }
      }, error = function(e) {
        cat("      Error:", conditionMessage(e), "\n")
        return(NULL)
      })
      
      if (is.null(model_fit)) {
        cat("      ", model, " model training failed. Skipping.\n")
        flush.console()
        next
      }
      
      model_pred <- tryCatch({
        predict(model_fit, newdata = test_data)
      }, error = function(e) {
        cat("      Prediction error:", conditionMessage(e), "\n")
        return(NULL)
      })
      
      if (is.null(model_pred)) {
        cat("      ", model, " model prediction failed. Skipping.\n")
        flush.console()
        next
      }
      
      model_actual <- test_data[[target]]
      
      # Calculate and store performance metrics
      Model_R2[model, target] <- tryCatch(cor(model_pred, model_actual)^2, error = function(e) NA)
      Model_RMSE[model, target] <- tryCatch(RMSE(model_pred, model_actual), error = function(e) NA)
      Model_MSE[model, target] <- tryCatch(mean((model_pred - model_actual)^2), error = function(e) NA)
      
      # Calculate and store variable importance (IMP)
      cat("      Calculating variable importance (IMP) for", model, "model...\n")
      flush.console()
      
      var_imp <- NULL
      if (model == "RF") {
          var_imp <- tryCatch({
            model_fit$finalModel$importance[env_vars, 'IncNodePurity']
          }, error = function(e) NULL)
      } else { # For GBM, SVR, KNN, GLM, NN
          var_imp <- compute_feature_imp_iml(model_fit, train_data, env_vars, target)
      }
      
      if (!is.null(var_imp)) {
        feature_importance_IMP[[model]][, target] <- var_imp
        cat("      IMP calculation for", model, "model successful.\n")
      } else {
        cat("      IMP calculation for", model, "model failed. Skipping save.\n")
      }

      # Calculate and store SHAP values
      cat("      Calculating SHAP values for", model, "model...\n")
      flush.console()
      shap_imp <- compute_shap(model_fit, train_data, env_vars, target, nsim = 100)
      
      if (!is.null(shap_imp)) {
        feature_importance_SHAP[[model]][, target] <- shap_imp
        cat("      SHAP value calculation for", model, "model successful.\n")
      } else {
        cat("      SHAP value calculation for", model, "model failed. Assigning NA.\n")
        feature_importance_SHAP[[model]][, target] <- NA
      }
      
      cat("    Finished training", model, "model for target variable:", target, ".\n\n")
      flush.console()
    }
  }
  
  # Define the save path for the current seed
  current_save_path <- file.path(output_root_path, seed_str, "/")
  
  # Create the save directory if it doesn't exist
  if (!dir.exists(current_save_path)) {
    dir.create(current_save_path, recursive = TRUE)
  }
  
  # Save all results for the current seed
  save_results(feature_importance_IMP, feature_importance_SHAP, Model_R2, Model_RMSE, Model_MSE, current_save_path, models)
  
  # Print a summary of the metrics used
  cat("  Summary of feature importance metrics used for this seed:\n")
  for (model in models) {
    cat("  Model:", model, "\n")
    if (model == "RF") {
      cat("    IMP - Variable importance calculated using IncNodePurity from the 'randomForest' package.\n")
    } else {
      cat("    IMP - Permutation importance calculated using FeatureImp from the 'iml' package.\n")
    }
    cat("    SHAP - SHAP values calculated using the 'fastshap' package.\n")
  }
  
  cat("=== Finished Seed:", seed_str, "===\n\n")
}

cat("Processing for all seeds is complete!\n")
