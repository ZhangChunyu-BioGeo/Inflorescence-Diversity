# =========================================================================
# Visualize Results of Multiple Machine Learning Models
# Target Variable: H_Rank
# Last updated: 2025-10-14
# =========================================================================

# Set the primary working directory for the project.
setwd("D:/R_work/Inflorescence")

# Define the target variable for this visualization script
VAR_T <- 'H_Rank'

# Load required packages
library(reshape2)
library(ggplot2)
library(dplyr)
library(cowplot)

# Define the list of models
models <- c("RF", "SVR", "KNN", "NN", "GLM", "GBM")

# Define the list of feature importance metrics
IMP_list <- c('IMP', 'SHAP')

# Create a list to store the plot objects
plot_list <- list()

# Loop through each importance metric (IMP and SHAP)
for (IMP in IMP_list) {
  
  # Initialize lists to store feature importance values and p-values for each model
  imp_list <- list()
  pval_list <- list()
  
  # Read the feature importance and p-value data for each model
  for (model in models) {
    # Construct file paths
    imp_file_path <- paste0("MLModel/Result/", model, "_", IMP, ".csv")
    pval_file_path <- paste0("MLModel/Result/", model, "_", IMP, "_p.csv")
    
    # Check if files exist
    if (!file.exists(imp_file_path)) {
      stop("File does not exist: ", imp_file_path)
    }
    if (!file.exists(pval_file_path)) {
      stop("File does not exist: ", pval_file_path)
    }
    
    # Read importance CSV file
    imp_df <- read.csv(imp_file_path, row.names = 1, check.names = FALSE)
    
    # Ensure the target variable column exists
    if (!(VAR_T %in% colnames(imp_df))) {
      stop("Target variable '", VAR_T, "' not found in file: ", imp_file_path)
    }
    
    # Extract importance values for the target variable
    imp_values <- imp_df[, VAR_T, drop = FALSE]
    imp_values <- as.numeric(as.character(imp_values[,1]))
    names(imp_values) <- rownames(imp_df)
    
    # Normalize importance values: absolute value then min-max scaling to [0, 1]
    imp_values_abs <- abs(imp_values)
    max_value <- max(imp_values_abs, na.rm = TRUE)
    min_value <- min(imp_values_abs, na.rm = TRUE)
    
    if (max_value == min_value) {
      # If all values are the same, normalize to 0.5 to avoid division by zero
      imp_values_norm <- rep(0.5, length(imp_values_abs))
    } else {
      imp_values_norm <- (imp_values_abs - min_value) / (max_value - min_value)
    }
    
    imp_list[[model]] <- imp_values_norm
    
    # Read the corresponding p-value CSV file
    pval_df <- read.csv(pval_file_path, row.names = 1, check.names = FALSE)
    
    if (!(VAR_T %in% colnames(pval_df))) {
      stop("Target variable '", VAR_T, "' not found in p-value file: ", pval_file_path)
    }
    
    # Extract p-values for the target variable
    pval_values <- pval_df[, VAR_T, drop = FALSE]
    pval_values <- as.numeric(as.character(pval_values[,1]))
    names(pval_values) <- rownames(pval_df)
    
    pval_list[[model]] <- pval_values
  }
  
  # Get predictor variable names from the RF model as a reference
  predictor_vars <- names(imp_list[['RF']])
  
  # Sort predictors based on the normalized importance values from the RF model
  rf_imp_values <- imp_list[['RF']]
  sorted_predictors <- names(sort(rf_imp_values, decreasing = TRUE))
  
  # Create a matrix to hold importance values for all models
  imp_matrix <- data.frame(matrix(nrow = length(models), ncol = length(sorted_predictors)))
  colnames(imp_matrix) <- sorted_predictors
  rownames(imp_matrix) <- models
  
  for (model in models) {
    imp_matrix[model, ] <- imp_list[[model]][sorted_predictors]
  }
  
  # Convert the matrix to a long format for ggplot2
  imp_long <- reshape2::melt(as.matrix(imp_matrix))
  colnames(imp_long) <- c("Model", "Variable", "IMP_Value")
  imp_long$IMP_Value <- as.numeric(as.character(imp_long$IMP_Value))
  
  # Set factor levels to ensure correct ordering in the plot
  imp_long$Model <- factor(imp_long$Model, levels = rev(models))
  imp_long$Variable <- factor(imp_long$Variable, levels = sorted_predictors)
  
  # Calculate within-model ranks for each variable
  imp_long <- imp_long %>%
    group_by(Model) %>%
    arrange(desc(IMP_Value)) %>%
    mutate(Rank = row_number()) %>%
    ungroup()
  
  # Create a long dataframe for p-values
  pval_long <- data.frame(
    Model = rep(models, each = length(sorted_predictors)),
    Variable = rep(sorted_predictors, times = length(models)),
    p_value = NA_real_
  )
  for (model in models) {
    pval_long$p_value[pval_long$Model == model] <- pval_list[[model]][sorted_predictors]
  }
  
  # Merge p-values into the main importance dataframe
  imp_long <- left_join(imp_long, pval_long, by = c("Model", "Variable"))
  
  # Define custom significance symbols and sizes based on p-value and rank
  imp_long <- imp_long %>%
    mutate(
      Top2 = Rank <= 2,
      Significance_Symbol = case_when(
        Top2 & p_value < 0.001 ~ "\u25CF",            # Top 2 and highly significant: Solid circle
        Top2 & p_value >= 0.001 ~ "\u25CB",           # Top 2 but less significant: Hollow circle
        !Top2 & p_value < 0.001 ~ "\u25CF",           # Not Top 2, but highly significant: Solid circle
        !Top2 & p_value < 0.01 ~ "\u25CB",            # Not Top 2, p < 0.01: Hollow circle
        !Top2 & p_value < 0.05 ~ "*",                 # Not Top 2, p < 0.05: Asterisk
        TRUE ~ "ns"                                   # Not significant
      ),
      Circle_Size = case_when(
        Significance_Symbol %in% c("\u25CF", "\u25CB") & Rank == 1 ~ 5.0, # Rank 1 circle size
        Significance_Symbol %in% c("\u25CF", "\u25CB") & Rank == 2 ~ 3.5, # Rank 2 circle size
        Significance_Symbol %in% c("\u25CF", "\u25CB") ~ 1.5,             # Other circles
        TRUE ~ NA_real_                               # No size for text symbols
      )
    )
  
  # --- Create the Heatmap Plot ---
  p <- ggplot(imp_long, aes(x = Variable, y = Model, fill = IMP_Value)) +
    geom_tile(color = "black", size = 0.2) +
    scale_fill_gradient(low = "white", high = "#F3A261") +
    # Add circles for high significance/rank
    geom_point(data = . %>% filter(!is.na(Circle_Size)),
               aes(size = Circle_Size, shape = Significance_Symbol),
               color = "#040676") +
    scale_size_identity() +
    scale_shape_manual(values = c("\u25CF" = 16, "\u25CB" = 1)) + # Map symbols to solid/hollow points
    # Add text for lower significance
    geom_text(data = . %>% filter(Significance_Symbol %in% c("*", "ns")),
              aes(label = Significance_Symbol),
              size = 4, color = "#040676") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = -45, hjust = 0),
      plot.margin = unit(c(1, 0, 1, 0), "cm"),
      legend.position = "none",
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    scale_y_discrete(limits = rev(models), position = "right") +
    scale_x_discrete(limits = sorted_predictors) +
    labs(x = "", y = "", fill = IMP)
  
  # Conditional formatting: only show y-axis labels for the SHAP plot
  if (IMP == "SHAP") {
    p <- p + theme(axis.text.y = element_text(face = "bold"))
  } else {
    p <- p + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  }
  
  # Store the finished plot in the list
  plot_list[[IMP]] <- p
}

# --- Create the Model Performance Plot ---

# Load model performance data (R2 and RMSE) and their p-values
Model_R2 <- read.csv("MLModel/Result/Model_R2.csv", row.names = 1, check.names = FALSE)
Model_RMSE <- read.csv("MLModel/Result/Model_RMSE.csv", row.names = 1, check.names = FALSE)
Model_R2_p <- read.csv("MLModel/Result/Model_R2_p.csv", row.names = 1, check.names = FALSE)
Model_RMSE_p <- read.csv("MLModel/Result/Model_RMSE_p.csv", row.names = 1, check.names = FALSE)

# Create a dataframe for plotting performance
performance_df <- data.frame(
  Model = factor(models, levels = rev(models)),
  R2 = Model_R2[models, VAR_T],
  RMSE = Model_RMSE[models, VAR_T],
  R2_p = Model_R2_p[models, VAR_T],
  RMSE_p = Model_RMSE_p[models, VAR_T]
)

# Convert performance values to numeric
performance_df <- performance_df %>%
  mutate(across(c(R2, RMSE, R2_p, RMSE_p), ~as.numeric(as.character(.))))

# Add significance labels for R2
performance_df <- performance_df %>%
  mutate(
    R2_Significance = case_when(
      R2_p < 0.001 ~ "***", R2_p < 0.01 ~ "**", R2_p < 0.05 ~ "*", TRUE ~ "ns"
    )
  )

# Scale RMSE to fit on the 0-1 range for plotting
performance_df$RMSE_Scaled <- pmin(performance_df$RMSE, 1.0) 

# Create the horizontal bar and point chart for performance
p3 <- ggplot(performance_df, aes(y = Model)) +
  # R2 bars
  geom_bar(aes(x = R2, fill = R2), stat = "identity", width = 0.6, color = "black", size = 0.4) +
  scale_fill_gradientn(colors = c("#FFF6D9", "#9BCCE3", "#264653"), limits = c(0.25, 0.75)) +
  # RMSE points and segments (plotted on a reversed scale)
  geom_point(aes(x = 1 - RMSE_Scaled), color = "black", size = 2) +
  geom_segment(aes(x = 1 - RMSE_Scaled, xend = 1 - RMSE_Scaled, 
                   y = as.numeric(Model) - 0.3, yend = as.numeric(Model) + 0.3), 
               color = "black", size = 0.4) +
  # Add R2 significance labels
  geom_text(aes(x = R2 + 0.02, label = R2_Significance), color = "#040676", size = 5, hjust = 0) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = unit(c(1, 1, 1, 0), "cm"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 0.7),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  # Define primary (R2) and secondary (RMSE) x-axes
  scale_x_continuous(
    name = "R²",
    sec.axis = sec_axis(~(1 - .), name = "RMSE") 
  )

# --- Combine All Plots ---
# Use cowplot to arrange the two importance heatmaps and the performance chart
combined_plot <- plot_grid(
  plot_list[['IMP']],
  plot_list[['SHAP']],
  p3,
  ncol = 3,
  align = 'h',
  rel_widths = c(2, 2, 1) # Adjust relative widths of the plots
)

# Display the final combined plot
print(combined_plot)
