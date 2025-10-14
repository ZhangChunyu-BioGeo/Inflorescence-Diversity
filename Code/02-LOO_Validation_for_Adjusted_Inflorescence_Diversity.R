# =========================================================================
# Leave-One-Out (LOO) Cross-Validation and Visualization for H_Rank
# Last updated: 2025-10-14
# =========================================================================

# Set the primary working directory for the project.
setwd("D:/R_work/Inflorescence")


# =========================================================================
# Part 1: Leave-One-Out (LOO) Calculation for Inflorescence Diversity (H_Rank)
# =========================================================================

# 1. Load Packages
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# 2. Define File Paths and Constants
# --- Input file ---
input_file <- "Geo_Data_Trait_02_Adjusted.csv"

# --- Intermediate output file for LOO results ---
loo_result_file <- "H_Rank_LOO/Geo_Data_H_Rank_LOO_Result.csv"

# Define column names for the nine proportion categories
proportion_cols <- c("So_Pi", "Df_Pi", "Co_Pi", "Sk_Pi", "Um_Pi", "Ra_Pi", "Cp_Pi", "Pa_Pi", "Ck_Pi")

# Visualization parameters
color_palette <- c("#56AA5D", "#f4f1f4", "#8D579C") 
category_order <- c("NoSo", "NoDf", "NoCo", "NoSk", "NoUm", "NoRa", "NoCp", "NoPa", "NoCk")
plot_category_order <- rev(category_order)

# 3. Define Function to Calculate H_Rank
calculate_h_rank <- function(df, prop_cols) {
  
  # Select the required proportion columns from the dataframe
  proportions <- df[, prop_cols, drop = FALSE]
  
  # Re-normalize proportions to ensure each row sums to 1
  row_sums <- rowSums(proportions)
  renormalized_props <- proportions / ifelse(row_sums == 0, 1, row_sums)
  
  # Calculate the Shannon diversity index (H)
  p_log_p <- renormalized_props * log(renormalized_props)
  p_log_p[is.na(p_log_p)] <- 0
  
  H <- -rowSums(p_log_p)
  
  # Get the total number of grid cells (G)
  G <- nrow(df)
  
  # Calculate the rank (R)
  R <- rank(H, ties.method = "average")
  
  # Normalize the rank to get H_Rank
  H_Rank <- R / (G + 1)
  
  return(H_Rank)
}

# 4. Perform LOO Calculation
tryCatch({
  data_for_loo <- read_csv(input_file)
}, error = function(e) {
  stop("Failed to read input file. Check if the path is correct and the file exists: ", input_file)
})

# Initialize a dataframe to store the results
loo_calc_results <- data.frame(Grid_ID = data_for_loo$Grid_ID)

# Calculate the baseline H_Rank using all nine categories
loo_calc_results$H_Rank <- calculate_h_rank(data_for_loo, proportion_cols)

# Loop through each category to perform the leave-one-out calculation
for (col_to_exclude in proportion_cols) {
  cols_to_use <- setdiff(proportion_cols, col_to_exclude)
  new_col_name <- paste0("H_Rank_No", sub("_Pi", "", col_to_exclude))
  cat(paste("  Calculating:", new_col_name, "...\n"))
  loo_calc_results[[new_col_name]] <- calculate_h_rank(data_for_loo, cols_to_use)
}

# 5. Save LOO Calculation Results
# Ensure the output directory exists before writing the file.
output_dir <- dirname(loo_result_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

tryCatch({
  write_csv(loo_calc_results, loo_result_file)
  cat(paste("LOO calculation completed successfully! Results saved to:\n", loo_result_file, "\n"))
}, error = function(e) {
  stop("Failed to write output file. Check if the path exists and write permissions are granted.")
})

cat("\nPreview of LOO results:\n")
print(head(loo_calc_results))


# =========================================================================
# Part 2: Visualization of Leave-One-Out (LOO) Results
# =========================================================================

# 6. Load Data for Visualization
tryCatch({
  vis_loo_results <- read_csv(loo_result_file)
}, error = function(e) { stop("Failed to read LOO results file: ", loo_result_file) })

tryCatch({
  vis_original_data <- read_csv(input_file, col_select = c("Grid_ID", "Species_N_Inf"))
}, error = function(e) { stop("Failed to read original data file: ", input_file) })

vis_loo_results <- left_join(vis_loo_results, vis_original_data, by = "Grid_ID")

# Create 5 quantile groups based on species number (Species_N_Inf)
quantile_breaks <- unique(quantile(vis_loo_results$Species_N_Inf, probs = seq(0, 1, 0.2), na.rm = TRUE))
if (length(quantile_breaks) < 2) {
  vis_loo_results <- vis_loo_results %>% mutate(SR_Quantile = as.factor(Species_N_Inf))
} else {
  vis_loo_results <- vis_loo_results %>%
    mutate(SR_Quantile = cut(Species_N_Inf, breaks = quantile_breaks, include.lowest = TRUE))
}

# 7. Data Transformation and Metric Calculation
vis_results_long <- vis_loo_results %>%
  pivot_longer(
    cols = starts_with("H_Rank_No"), names_to = "LOO_Method_Raw", values_to = "LOO_H_Rank"
  ) %>%
  mutate(
    Difference = H_Rank - LOO_H_Rank,
    LOO_Method = sub("H_Rank_", "", LOO_Method_Raw)
  ) %>%
  mutate(LOO_Method = factor(LOO_Method, levels = plot_category_order))

# Calculate R-squared for each LOO scenario
vis_scatter_metrics <- vis_results_long %>%
  group_by(LOO_Method) %>%
  summarise(r_squared = cor(H_Rank, LOO_H_Rank)^2, .groups = 'drop') %>%
  mutate(r_squared_label = sprintf("R² = %.4f", r_squared))

# Calculate mean and standard deviation of the difference
vis_boxplot_metrics <- vis_results_long %>%
  group_by(LOO_Method) %>%
  summarise(mean_diff = mean(Difference), sd_diff = sd(Difference), .groups = 'drop') %>%
  mutate(metrics_label = sprintf("Mean = %+.4f\nSD = %.4f", mean_diff, sd_diff))

# Join the standard deviation back to the main long-format dataframe
vis_results_long <- left_join(vis_results_long, 
                               vis_boxplot_metrics %>% select(LOO_Method, sd_diff), 
                               by = "LOO_Method")

# 8. Visualization 1: Scatter Plot Matrix
scatter_5_colors <- c("#095624", "#9CCBA1", "#D1D1D1", "#BA9AC3", "#40004b")

vis_scatter_plot <- ggplot(vis_results_long, aes(x = H_Rank, y = LOO_H_Rank)) +
  geom_point(aes(color = SR_Quantile), alpha = 0.5, size = 1) +
  scale_color_manual(values = scatter_5_colors, name = "Species Number", na.translate = FALSE) +
  geom_abline(intercept = 0, slope = 1, color = "#D73831", linetype = "dashed", linewidth = 0.8) +
  geom_text(
    data = vis_scatter_metrics, aes(label = r_squared_label, y = LOO_Method),
    x = 0.05, y = 0.95, hjust = 0, vjust = 1, size = 3.5, fontface = "bold"
  ) +
  facet_wrap(~ factor(LOO_Method, levels = category_order), ncol = 3) +
  labs(
    x = "Original H_Rank (All 9 Categories)",
    y = "Leave-One-Out H_Rank (8 Categories)"
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5))) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_text(face="bold")
  )

# 9. Visualization 2: Box Plot of Differences
vis_outlier_data <- vis_results_long %>%
  group_by(LOO_Method) %>%
  mutate(
    q1 = quantile(Difference, 0.25, na.rm = TRUE),
    q3 = quantile(Difference, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    is_outlier = Difference < (q1 - 1.5 * iqr) | Difference > (q3 + 1.5 * iqr)
  ) %>%
  filter(is_outlier)

x_range <- range(vis_results_long$Difference, na.rm = TRUE)
text_x_pos <- x_range[2] + (x_range[2] - x_range[1]) * 0.05

vis_boxplot <- ggplot(vis_results_long, aes(y = LOO_Method, x = Difference)) +
  geom_boxplot(aes(fill = sd_diff), 
               color = "#262626", 
               linewidth = 0.3,
               outlier.shape = NA) +
  geom_point(data = vis_outlier_data, 
             aes(color = sd_diff), 
             alpha = 0.6, 
             shape = 19,
             size = 1.5) +
  scale_fill_gradientn(colors = color_palette, name = "Std. Deviation (SD)\n of Difference") +
  scale_color_gradientn(colors = color_palette, guide = "none") +
  geom_vline(xintercept = 0, color = "#262626", linetype = "dashed", linewidth = 0.8) +
  geom_text(
    data = vis_boxplot_metrics, aes(label = metrics_label),
    x = text_x_pos, hjust = 0,
    size = 3.5, lineheight = .9, fontface = "bold"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  labs(
    y = "Excluded Category",
    x = "Change in H_Rank"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 11, face = "bold"),
    legend.position = "right",
    legend.title = element_text(face="bold")
  )

# 10. Display Plots
print(vis_scatter_plot)
dev.new()
print(vis_boxplot)
