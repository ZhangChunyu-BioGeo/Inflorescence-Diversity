# =========================================================================
# Visualization of SCM (Stochastic Character Mapping) Results
# Last updated: 2025-10-14
# =========================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Set the primary working directory for the project.
setwd("D:/R_work/Inflorescence")

# -------------------------------------------------------------------------
# Part 1: Parameters and Setup
# -------------------------------------------------------------------------
scm_summary_file <- "Evolution_Model/SCM_Summary_HR2_ARD.csv"
climate_file <- "Evolution_Model/Aggregated_Climate_Data.csv"
output_dir <- "Evolution_Model/"

# --- Optional: Number of oldest bins to drop from visualizations ---
# This can help focus on more recent periods with potentially better data.
DROP_OLDEST_BINS_P1 <- 0 # For Figure 1 (Proportions + GMST)
DROP_OLDEST_BINS_P2 <- 3 # For Figure 2 (Transitions + Rates)
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Part 2: Data Loading and Preparation
# -------------------------------------------------------------------------
# Load SCM summary data
scm_data <- read.csv(scm_summary_file, stringsAsFactors = FALSE)
cat("SCM summary data loaded:", nrow(scm_data), "rows.\n")

# Load aggregated climate data
climate_data <- read.csv(climate_file, stringsAsFactors = FALSE)
cat("Aggregated climate data loaded:", nrow(climate_data), "rows.\n")

# Merge SCM and climate data by bin_id
combined_data <- merge(scm_data, 
                       climate_data[, c("bin_id", "GMST_05", "GMST_16", "GMST_50", "GMST_84", "GMST_95", "coverage_percent")], 
                       by = "bin_id", all.x = TRUE)

# Calculate 95% confidence intervals (mean ± 1.96 * SD)
combined_data <- combined_data %>%
  mutate(
    # 95% CI for proportions, bounded between 0 and 1
    prop_D_ci_lower = pmax(0, prop_D_mean - 1.96 * prop_D_sd),
    prop_D_ci_upper = pmin(1, prop_D_mean + 1.96 * prop_D_sd),
    prop_I_ci_lower = pmax(0, prop_I_mean - 1.96 * prop_I_sd),
    prop_I_ci_upper = pmin(1, prop_I_mean + 1.96 * prop_I_sd),
    
    # 95% CI for transition counts, bounded at a minimum of 0
    trans_DtoI_ci_lower = pmax(0, transitions_DtoI_mean - 1.96 * transitions_DtoI_sd),
    trans_DtoI_ci_upper = transitions_DtoI_mean + 1.96 * transitions_DtoI_sd,
    trans_ItoD_ci_lower = pmax(0, transitions_ItoD_mean - 1.96 * transitions_ItoD_sd),
    trans_ItoD_ci_upper = transitions_ItoD_mean + 1.96 * transitions_ItoD_sd,
    
    # 95% CI for transition rates, bounded between 0 and a plausible max (0.035)
    rate_DtoI_ci_lower = pmax(0, rate_DtoI_per_My_mean - 1.96 * rate_DtoI_per_My_sd),
    rate_DtoI_ci_upper = pmin(0.035, rate_DtoI_per_My_mean + 1.96 * rate_DtoI_per_My_sd),
    rate_ItoD_ci_lower = pmax(0, rate_ItoD_per_My_mean - 1.96 * rate_ItoD_per_My_sd),
    rate_ItoD_ci_upper = pmin(0.035, rate_ItoD_per_My_mean + 1.96 * rate_ItoD_per_My_sd)
  ) %>%
  arrange(midpoint_age)

# --- Define consistent color schemes ---
trait_colors <- c("D" = "#56AA5E", "I" = "#8B559B")
direction_colors <- c("D→I" = "#8B559B", "I→D" = "#56AA5E")


# -------------------------------------------------------------------------
# Part 3: Figure 1 - Trait Proportions and Paleoclimate (GMST)
# -------------------------------------------------------------------------

# Prepare data for Figure 1, dropping oldest bins if specified
p1_data <- combined_data
if (DROP_OLDEST_BINS_P1 > 0 && nrow(p1_data) > DROP_OLDEST_BINS_P1) {
  cat("--- Preparing data for Figure 1: Dropping", DROP_OLDEST_BINS_P1, "oldest bins ---\n")
  p1_data <- p1_data %>%
    arrange(desc(midpoint_age)) %>% # Sort by age, oldest first
    slice(-(1:DROP_OLDEST_BINS_P1)) %>% # Remove the top N rows
    arrange(midpoint_age) # Revert to original order
} else {
  cat("--- Preparing data for Figure 1: Using all", nrow(p1_data), "bins ---\n")
}

# --- Prepare GMST data for plotting ---
gmst_data <- p1_data %>%
  filter(!is.na(GMST_50)) %>%
  arrange(midpoint_age) %>%
  mutate(
    # Normalize GMST to a 0-1 scale based on a 0-45°C range for plotting
    gmst_05_norm = pmax(0, pmin(1, GMST_05 / 45)),
    gmst_16_norm = pmax(0, pmin(1, GMST_16 / 45)),
    gmst_50_norm = pmax(0, pmin(1, GMST_50 / 45)),
    gmst_84_norm = pmax(0, pmin(1, GMST_84 / 45)),
    gmst_95_norm = pmax(0, pmin(1, GMST_95 / 45))
  )

# Create dataframes for segments and polygons to give a continuous look
if(nrow(gmst_data) > 1) {
  # Median line (GMST_50)
  line_segments <- bind_rows(lapply(1:(nrow(gmst_data)-1), function(i) {
    data.frame(x = gmst_data$midpoint_age[i], y = gmst_data$gmst_50_norm[i],
               xend = gmst_data$midpoint_age[i+1], yend = gmst_data$gmst_50_norm[i+1],
               gmst_avg = (gmst_data$GMST_50[i] + gmst_data$GMST_50[i+1]) / 2)
  }))
  # Outer confidence interval (5th-95th percentile)
  outer_polygons <- bind_rows(lapply(1:(nrow(gmst_data)-1), function(i) {
    data.frame(x = c(gmst_data$midpoint_age[i], gmst_data$midpoint_age[i+1], gmst_data$midpoint_age[i+1], gmst_data$midpoint_age[i]),
               y = c(gmst_data$gmst_05_norm[i], gmst_data$gmst_05_norm[i+1], gmst_data$gmst_95_norm[i+1], gmst_data$gmst_95_norm[i]),
               group = i, gmst_avg = (gmst_data$GMST_50[i] + gmst_data$GMST_50[i+1]) / 2)
  }))
  # Inner confidence interval (16th-84th percentile)
  inner_polygons <- bind_rows(lapply(1:(nrow(gmst_data)-1), function(i) {
    data.frame(x = c(gmst_data$midpoint_age[i], gmst_data$midpoint_age[i+1], gmst_data$midpoint_age[i+1], gmst_data$midpoint_age[i]),
               y = c(gmst_data$gmst_16_norm[i], gmst_data$gmst_16_norm[i+1], gmst_data$gmst_84_norm[i+1], gmst_data$gmst_84_norm[i]),
               group = i, gmst_avg = (gmst_data$GMST_50[i] + gmst_data$GMST_50[i+1]) / 2)
  }))
}

# --- Prepare data for stacked proportions ---
stacked_data <- p1_data %>%
  mutate(
    prop_I_bottom = 0,
    prop_I_top = prop_I_mean,
    prop_D_bottom = prop_I_mean,
    prop_D_top = 1,
    error_bar_y = prop_I_mean,
    error_bar_lower = prop_I_ci_lower,
    error_bar_upper = prop_I_ci_upper
  )

# Calculate GMST range for the color legend midpoint
gmst_mid <- mean(range(p1_data$GMST_50, na.rm = TRUE))

# --- Build Figure 1 ---
p1 <- ggplot(stacked_data, aes(x = midpoint_age)) +
  # Layer 1: Stacked ribbons for trait proportions
  geom_ribbon(aes(ymin = prop_I_bottom, ymax = prop_I_top), fill = trait_colors["I"], alpha = 0.5) +
  geom_ribbon(aes(ymin = prop_D_bottom, ymax = prop_D_top), fill = trait_colors["D"], alpha = 0.5) +
  # Layer 2: Error bars on the proportion boundary
  geom_errorbar(aes(ymin = error_bar_lower, ymax = error_bar_upper), width = 0, color = "#262626", size = 0.5, alpha = 0.95) +
  # Layer 3: GMST confidence intervals (polygons) and median (line)
  {if(exists("outer_polygons")) geom_polygon(data = outer_polygons, aes(x=x, y=y, group=group, fill=gmst_avg), alpha = 0.4)} +
  {if(exists("inner_polygons")) geom_polygon(data = inner_polygons, aes(x=x, y=y, group=group, fill=gmst_avg), alpha = 0.6)} +
  {if(exists("line_segments")) geom_segment(data = line_segments, aes(x=x, y=y, xend=xend, yend=yend, color=gmst_avg), size=1.2)} +
  # --- Scales and Axes ---
  scale_fill_gradient2(low = "#93B6D5", mid = "white", high = "#F0869A", midpoint = gmst_mid, name = "GMST (°C)", guide = guide_colorbar(title.position="top", title.hjust=0.5, barwidth=1, barheight=8)) +
  scale_color_gradient2(low = "#93B6D5", mid = "white", high = "#F0869A", midpoint = gmst_mid, guide = "none") +
  scale_x_reverse(name = "Evolutionary Time (Ma)", breaks = pretty_breaks(n = 8)) +
  scale_y_continuous(
    name = "Inflorescence Type Proportion", limits = c(0, 1), breaks = seq(0, 1, 0.2),
    sec.axis = sec_axis(trans = ~ . * 45, name = "Global Mean Surface Temperature (°C)", breaks = seq(0, 45, 5))
  ) +
  # --- Theme ---
  theme_classic() +
  theme(
    panel.border = element_rect(color="black", fill=NA, size=0.8), axis.line=element_line(),
    axis.text = element_text(size=11, color="black"), axis.title = element_text(size=12, color="black"),
    legend.position="right", legend.title=element_text(size=11), legend.text=element_text(size=10),
    legend.background = element_rect(color="black", size=0.5), legend.margin = margin(5,5,5,5)
  )

# -------------------------------------------------------------------------
# Part 4: Figure 2 - Transition Counts and Rates
# -------------------------------------------------------------------------

# Prepare data for Figure 2, dropping oldest bins if specified
p2_data <- combined_data
if (DROP_OLDEST_BINS_P2 > 0 && nrow(p2_data) > DROP_OLDEST_BINS_P2) {
  cat("--- Preparing data for Figure 2: Dropping", DROP_OLDEST_BINS_P2, "oldest bins ---\n")
  p2_data <- p2_data %>%
    arrange(desc(midpoint_age)) %>% # Sort by age, oldest first
    slice(-(1:DROP_OLDEST_BINS_P2)) %>% # Remove the top N rows
    arrange(midpoint_age) # Revert to original order
} else {
  cat("--- Preparing data for Figure 2: Using all", nrow(p2_data), "bins ---\n")
}

# --- Prepare data for plotting by pivoting to long format ---
# Transition counts
transition_counts_data <- p2_data %>%
  select(midpoint_age, starts_with("transitions_"), starts_with("trans_")) %>%
  pivot_longer(cols = c(transitions_DtoI_mean, transitions_ItoD_mean), names_to = "direction", values_to = "count") %>%
  mutate(
    direction_label = if_else(grepl("DtoI", direction), "D→I", "I→D"),
    ci_lower = if_else(grepl("DtoI", direction), trans_DtoI_ci_lower, trans_ItoD_ci_lower),
    ci_upper = if_else(grepl("DtoI", direction), trans_DtoI_ci_upper, trans_ItoD_ci_upper)
  )

# Transition rates
transition_rates_data <- p2_data %>%
  select(midpoint_age, starts_with("rate_")) %>%
  pivot_longer(cols = c(rate_DtoI_per_My_mean, rate_ItoD_per_My_mean), names_to = "direction", values_to = "rate") %>%
  mutate(
    direction_label = if_else(grepl("DtoI", direction), "D→I", "I→D"),
    ci_lower = if_else(grepl("DtoI", direction), rate_DtoI_ci_lower, rate_ItoD_ci_lower),
    ci_upper = if_else(grepl("DtoI", direction), rate_DtoI_ci_upper, rate_ItoD_ci_upper)
  )

# Determine scaling factor for the secondary axis
max_count <- max(transition_counts_data$ci_upper, na.rm = TRUE)
conversion_factor <- max_count / 0.035 # Map max count to a rate of 0.035
left_axis_max <- max_count * 1.05

# --- Build Figure 2 ---
p2 <- ggplot() +
  # Layer 1: Columns for transition counts
  geom_col(data = transition_counts_data, aes(x = midpoint_age, y = count, fill = direction_label),
           alpha = 0.5, width = 3.5, position = position_dodge(width = 3.5)) +
  # Layer 2: Error bars for counts
  geom_errorbar(data = transition_counts_data, aes(x = midpoint_age, ymin = ci_lower, ymax = ci_upper, color = direction_label),
                width = 0, position = position_dodge(width = 3.5), size = 0.5) +
  # Layer 3: Lines, points, and ribbons for transition rates (scaled)
  geom_line(data = transition_rates_data, aes(x = midpoint_age, y = rate * conversion_factor, color = direction_label), size = 1) +
  geom_point(data = transition_rates_data, aes(x = midpoint_age, y = rate * conversion_factor, color = direction_label), size = 1.5) +
  geom_ribbon(data = transition_rates_data, aes(x = midpoint_age, ymin = ci_lower * conversion_factor, ymax = ci_upper * conversion_factor, fill = direction_label), alpha = 0.2) +
  # --- Scales and Axes ---
  scale_x_reverse(name = "Evolutionary Time (Ma)", breaks = pretty_breaks(n = 8)) +
  scale_y_continuous(
    name = "Number of Transitions", limits = c(0, left_axis_max), breaks = pretty_breaks(n = 6),
    sec.axis = sec_axis(trans = ~ . / conversion_factor, name = "Transition Rate (per Myr)", breaks = seq(0, 0.035, 0.005))
  ) +
  scale_fill_manual(values = direction_colors, name = "Transition Direction") +
  scale_color_manual(values = direction_colors, name = "Transition Direction") +
  # --- Theme ---
  theme_classic() +
  theme(
    panel.border = element_rect(color="black", fill=NA, size=0.8), axis.line=element_line(),
    axis.text = element_text(size=11, color="black"), axis.title = element_text(size=12, color="black"),
    legend.position = "top", legend.title = element_text(size = 11), legend.text = element_text(size = 10)
  )

# -------------------------------------------------------------------------
# Part 5: Display and Save Plots
# -------------------------------------------------------------------------

# Display the plots in the R session
print(p1)

dev.new()
print(p2)
