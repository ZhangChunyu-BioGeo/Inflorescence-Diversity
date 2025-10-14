# =========================================================================
# Regression Analysis of GMST on Transition Rates and Proportions
# Last updated: 2025-10-14
# =========================================================================

# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
# library(broom) # broom is not used in the original code, so it's commented out.

# Set the primary working directory for the project.
setwd("D:/R_work/Inflorescence")

# -------------------------------------------------------------------------
# Part 1: Data Loading and Merging
# -------------------------------------------------------------------------
scm_summary_file <- "Evolution_Model/SCM_Summary_HR2_ARD.csv"
climate_file <- "Evolution_Model/Aggregated_Climate_Data.csv"

# Read the datasets
scm_data <- read.csv(scm_summary_file, stringsAsFactors = FALSE)
climate_data <- read.csv(climate_file, stringsAsFactors = FALSE)

# Merge SCM and climate data
combined_data <- merge(scm_data, 
                       climate_data[, c("bin_id", "GMST_50")], 
                       by = "bin_id", all.x = TRUE)
combined_data <- combined_data %>% arrange(midpoint_age)


# -------------------------------------------------------------------------
# Part 2: Exclude the Oldest Bins from Analysis
# -------------------------------------------------------------------------
cat("--- Excluding the 3 oldest bins from the analysis ---\n")
excluded_bins_count <- 3
total_bins <- nrow(combined_data)

if (excluded_bins_count >= total_bins) {
  stop("Error: Number of bins to exclude is greater than or equal to the total number of bins.")
}

# Arrange by age (descending) and remove the first N rows
combined_data_filtered <- combined_data %>%
  arrange(desc(midpoint_age)) %>%
  slice(-(1:excluded_bins_count)) %>%
  arrange(midpoint_age) 

if (nrow(combined_data_filtered) > 0) {
  cat("Analysis will be performed on", nrow(combined_data_filtered), "bins.\n")
  cat("Age range after exclusion:", round(min(combined_data_filtered$midpoint_age), 1), 
      "to", round(max(combined_data_filtered$midpoint_age), 1), "Ma\n")
}


# -------------------------------------------------------------------------
# Part 3: Variable Selection and Data Preparation
# -------------------------------------------------------------------------
# Define the variables of interest and their labels for plotting
prop_var <- "prop_I_mean"
prop_label <- "Indeterminate Proportion"
rate_var <- "rate_DtoI_per_My_mean"
rate_label <- "D→I Transition Rate"

# Create the final analysis dataset, removing any rows with missing data
analysis_data <- combined_data_filtered %>%
  filter(!is.na(GMST_50) & !is.na(.data[[rate_var]]) & !is.na(.data[[prop_var]])) %>%
  select(bin_id, midpoint_age, GMST_50, all_of(c(rate_var, prop_var))) %>%
  rename(
    rate_value = !!sym(rate_var),
    prop_value = !!sym(prop_var)
  )

cat("Final number of data points for regression:", nrow(analysis_data), "\n")


# -------------------------------------------------------------------------
# Part 4: Linear Regression Analysis
# -------------------------------------------------------------------------
# Initialize variables to store results
rate_r2 <- 0; rate_p <- 1
prop_r2 <- 0; prop_p <- 1

# Ensure there are enough data points to run a regression
if (nrow(analysis_data) >= 4) {
  tryCatch({
    # --- Linear regression for Proportion vs. GMST ---
    lm_prop <- lm(prop_value ~ GMST_50, data = analysis_data)
    prop_summary <- summary(lm_prop)
    prop_r2 <- prop_summary$r.squared
    prop_p <- coef(prop_summary)["GMST_50", "Pr(>|t|)"] # p-value for the GMST coefficient
    cat("Proportion ~ GMST linear model - R²:", round(prop_r2, 4), "| p-value:", signif(prop_p, 4), "\n")

    # --- Linear regression for Rate vs. GMST ---
    lm_rate <- lm(rate_value ~ GMST_50, data = analysis_data)
    rate_summary <- summary(lm_rate)
    rate_r2 <- rate_summary$r.squared
    rate_p <- coef(rate_summary)["GMST_50", "Pr(>|t|)"] # p-value for the GMST coefficient
    cat("Rate ~ GMST linear model       - R²:", round(rate_r2, 4), "| p-value:", signif(rate_p, 4), "\n")

  }, error = function(e) {
    cat("An error occurred during linear regression:", e$message, "\n")
  })
} else {
  cat("Skipping regression due to insufficient data (fewer than 4 points).\n")
}


# -------------------------------------------------------------------------
# Part 5: Visualization with Dual Y-Axis
# -------------------------------------------------------------------------
# Calculate scaling parameters to map the rate axis to the proportion axis
rate_range <- range(analysis_data$rate_value, na.rm = TRUE)
prop_range <- range(analysis_data$prop_value, na.rm = TRUE)

# To avoid division by zero or extreme scaling if range is small
if (diff(rate_range) < 1e-6) rate_range[2] <- rate_range[1] + 1e-6
if (diff(prop_range) < 1e-6) prop_range[2] <- prop_range[1] + 1e-6

scale_factor <- diff(prop_range) / diff(rate_range)
offset <- prop_range[1] - rate_range[1] * scale_factor

# --- Create the plot ---
p_reg <- ggplot(analysis_data, aes(x = GMST_50)) +
  # Plot for Proportion
  geom_point(aes(y = prop_value, color = "proportion"), size = 3, alpha = 0.7, shape = 16) +
  geom_smooth(aes(y = prop_value, color = "proportion", fill = "proportion"), 
              method = "lm", se = TRUE, alpha = 0.2, linewidth = 1.2) +
  
  # Plot for Rate (scaled to the proportion axis)
  geom_point(aes(y = rate_value * scale_factor + offset, color = "rate"), size = 3, alpha = 0.7, shape = 17) + # triangle shape
  geom_smooth(aes(y = rate_value * scale_factor + offset, color = "rate", fill = "rate"), 
              method = "lm", se = TRUE, alpha = 0.2, linewidth = 1.2, linetype = "dashed") +
  
  # --- Scales and Labels ---
  scale_color_manual(name = "Variable",
                     values = c("rate" = "#686868", "proportion" = "#8B559B"),
                     labels = c("rate" = rate_label, "proportion" = prop_label)) +
  scale_fill_manual(name = "Variable",
                    values = c("rate" = "#686868", "proportion" = "#8B559B"),
                    labels = c("rate" = rate_label, "proportion" = prop_label)) +
  
  scale_y_continuous(
    name = prop_label,
    sec.axis = sec_axis(~ (. - offset) / scale_factor, name = rate_label)
  ) +
  
  labs(x = "Global Mean Surface Temperature (°C)") +
  
  # --- Theme ---
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.position = "bottom"
  )

# --- Add regression statistics as annotations on the plot ---
gmst_range <- range(analysis_data$GMST_50, na.rm = TRUE)
annotation_x <- gmst_range[1]

p_reg <- p_reg +
  # Annotation for Proportion
  annotate("text", x = annotation_x, y = prop_range[2], 
           label = paste0("R² = ", round(prop_r2, 3), ", p = ", signif(prop_p, 3)),
           hjust = -0.1, vjust = 1.5, color = "#8B559B", size = 4, fontface = "bold") +
  # Annotation for Rate
  annotate("text", x = annotation_x, y = prop_range[1], 
           label = paste0("R² = ", round(rate_r2, 3), ", p = ", signif(rate_p, 3)),
           hjust = -0.1, vjust = -0.5, color = "#686868", size = 4, fontface = "bold")

# --- Display the final plot ---
print(p_reg)
