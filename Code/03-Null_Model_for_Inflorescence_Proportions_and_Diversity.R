# =========================================================================
# Null Model Simulation and Visualization for Inflorescence Proportions and Diversity
# Last updated: 2025-10-14
# =========================================================================

# Load packages
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(broom)
library(patchwork)
library(factoextra)
library(cowplot)

# Set the primary working directory for the project.
setwd("D:/R_work/Inflorescence")

# =========================================================================
# Part 1: Null Model Simulation
# =========================================================================

# 1. Data Import
# Import the comprehensive dataset containing traits and environmental data.
df <- read.csv(file = "Geo_Data_Trait_and_Enveriment.csv", header = TRUE)

# 2. Calculate Global Trait Proportions
total_species <- sum(df$Species_N_Inf)
global_counts <- df %>%
  summarize(
    So = sum(Species_N_Inf * So_Pi, na.rm = TRUE),
    Df = sum(Species_N_Inf * Df_Pi, na.rm = TRUE),
    Co = sum(Species_N_Inf * Co_Pi, na.rm = TRUE),
    Sk = sum(Species_N_Inf * Sk_Pi, na.rm = TRUE),
    Um = sum(Species_N_Inf * Um_Pi, na.rm = TRUE),
    Ra = sum(Species_N_Inf * Ra_Pi, na.rm = TRUE),
    Cp = sum(Species_N_Inf * Cp_Pi, na.rm = TRUE),
    Pa = sum(Species_N_Inf * Pa_Pi, na.rm = TRUE),
    Ck = sum(Species_N_Inf * Ck_Pi, na.rm = TRUE)
  )
global_proportions <- global_counts %>%
  mutate(across(everything(), ~ . / total_species))

# Convert to a probability vector for the simulation
prob_vector <- as.numeric(global_proportions)
names(prob_vector) <- c("So", "Df", "Co", "Sk", "Um", "Ra", "Cp", "Pa", "Ck")
prob_vector <- prob_vector / sum(prob_vector)

# 3. Calculate Parameters for Bayesian Adjustment
total_Indeterminate_species <- sum(df$Indeterminate * df$Species_N_Inf, na.rm = TRUE)
total_species_obs <- sum(df$Species_N_Inf, na.rm = TRUE)
p0 <- total_Indeterminate_species / total_species_obs

mean_p <- mean(df$Indeterminate, na.rm = TRUE)
var_p <- var(df$Indeterminate, na.rm = TRUE)
alpha <- max(mean_p * (1 - mean_p) / var_p - 1, 0.01)

# 4. Set Simulation Parameters
num_sim <- 1000
num_sites <- nrow(df)
attributes <- c("So", "Df", "Co", "Sk", "Um", "Ra", "Cp", "Pa", "Ck")

# Initialize matrices to store simulation results
sim_indeterminate <- matrix(0, nrow = num_sites, ncol = num_sim)
sim_indeterminate_adj <- matrix(0, nrow = num_sites, ncol = num_sim)
sim_H <- matrix(0, nrow = num_sites, ncol = num_sim)
sim_H_Rank <- matrix(0, nrow = num_sites, ncol = num_sim)

# Define function to compute Shannon's diversity index
compute_entropy <- function(proportions) {
  proportions <- proportions[proportions > 0]
  if(length(proportions) == 0) return(0)
  -sum(proportions * log(proportions))
}

# 5. Run the Simulation
set.seed(123)
for(i in 1:num_sites){
  n_species <- df$Species_N_Inf[i]
  counts <- rmultinom(n = num_sim, size = n_species, prob = prob_vector)
  props <- counts / n_species
  
  # Calculate proportion of indeterminate inflorescences
  sim_indeterminate[i, ] <- colSums(props[c("Co", "Sk", "Um", "Ra", "Cp", "Pa", "Ck"), , drop = FALSE])
  
  # Calculate Bayesian-adjusted proportion of indeterminate inflorescences
  sim_indeterminate_adj[i, ] <- (sim_indeterminate[i, ] * df$Species_N_Inf[i] + p0 * alpha) / (df$Species_N_Inf[i] + alpha)
  
  # Calculate Shannon diversity
  sim_H[i, ] <- apply(props, 2, compute_entropy)
}

# 6. Rank-Transform Simulated H Values
for(j in 1:num_sim){
  ranks <- rank(sim_H[, j], ties.method = "average")
  sim_H_Rank[, j] <- ranks / (num_sites + 1)
}

# 7. Calculate Summary Statistics for Simulation Results
indeterminate_mean <- rowMeans(sim_indeterminate)
indeterminate_sd <- apply(sim_indeterminate, 1, sd)

indeterminate_adj_mean <- rowMeans(sim_indeterminate_adj)
indeterminate_adj_sd <- apply(sim_indeterminate_adj, 1, sd)

H_Rank_mean <- rowMeans(sim_H_Rank)
H_Rank_sd <- apply(sim_H_Rank, 1, sd)

determinate_mean <- 1 - indeterminate_adj_mean
determinate_sd <- indeterminate_adj_sd # Standard deviation remains the same

# 8. Append Simulation Results to the Dataframe
df <- df %>%
  mutate(
    Indeterminate_mean = indeterminate_mean,
    Indeterminate_sd = indeterminate_sd,
    Indeterminate_adj_mean = indeterminate_adj_mean,
    Indeterminate_adj_sd = indeterminate_adj_sd,
    Determinate_mean = determinate_mean,
    Determinate_sd = determinate_sd,
    H_Rank_mean = H_Rank_mean,
    H_Rank_sd = H_Rank_sd
  )

# 9. Save Simulation Results
output_dir <- "Null_Model_Sim"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
write.csv(df, file = "Null_Model_Sim/Geo_Data_Trait_Null_Model.csv", row.names = FALSE)


# =========================================================================
# Part 2: Visualization
# =========================================================================

# Define environmental data columns for convenience
env_data_subset <- df[, 22:46]

# --- Plot 1: H_Rank vs. PC1 (p1) ---
# Data preparation
df_p1_long <- df %>%
  select(PC1, H_Rank_mean, H_Rank) %>%
  rename(NullModel = H_Rank_mean, Observed = H_Rank) %>%
  pivot_longer(cols = c(NullModel, Observed), names_to = "Group", values_to = "Value") %>%
  mutate(Group = factor(Group, levels = c("NullModel", "Observed")))

# Correct for log transformation: calculate shift to ensure positive arguments
min_pc1 <- min(df$PC1, na.rm = TRUE)
shift_p1 <- ifelse(min_pc1 <= 0, abs(min_pc1) + 1, 0)

# Model fitting and label generation
fit_null_p1 <- lm(Value ~ log(PC1 + shift_p1), data = df_p1_long %>% filter(Group == "NullModel"))
fit_obs_p1 <- lm(Value ~ poly(PC1, 2, raw = TRUE), data = df_p1_long %>% filter(Group == "Observed"))

label_null_p1 <- paste0("y = ", round(coef(fit_null_p1)[1], 2), " + ", round(coef(fit_null_p1)[2], 2),
                      " ln(x + ", round(shift_p1, 2), "), R² = ", round(summary(fit_null_p1)$r.squared, 3), "***")
label_obs_p1 <- paste0("y = ", round(coef(fit_obs_p1)[3], 2), "(x - ", round(-coef(fit_obs_p1)[2] / (2*coef(fit_obs_p1)[3]), 2),
                     ")² + ", round(coef(fit_obs_p1)[1] - coef(fit_obs_p1)[2]^2 / (4*coef(fit_obs_p1)[3]), 2),
                     ", R² = ", round(summary(fit_obs_p1)$r.squared, 3), "***")

# Plotting
p1 <- ggplot(data = df_p1_long, aes(x = PC1, y = Value, color = Group, fill = Group)) +
  geom_point(alpha = 0.1, size = 1.5) +
  geom_smooth(data = . %>% filter(Group == "NullModel"), method = "lm", formula = y ~ log(x + shift_p1), se = TRUE, alpha = 0.3) +
  geom_smooth(data = . %>% filter(Group == "Observed"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = TRUE, alpha = 0.3) +
  annotate("text", x = min(df$PC1, na.rm=TRUE), y = max(df_p1_long$Value, na.rm=TRUE), label = label_null_p1, hjust = 0, vjust = 1, size = 4, color = "black") +
  annotate("text", x = min(df$PC1, na.rm=TRUE), y = max(df_p1_long$Value, na.rm=TRUE) * 0.92, label = label_obs_p1, hjust = 0, vjust = 1, size = 4, color = "black") +
  scale_color_manual(values = c("NullModel" = "#707070", "Observed" = "#8B559B")) +
  scale_fill_manual(values = c("NullModel" = "#707070", "Observed" = "#8B559B")) +
  labs(x = "PC1", y = "H_Rank") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = c(0.2, 0.75),
    legend.title = element_blank(),
    legend.background = element_blank()
  ) +
  # Set Y-axis limits to match the data range
  coord_cartesian(ylim = range(df_p1_long$Value, na.rm = TRUE))

# --- Plot 2: Indeterminate_adj vs. PC1 (p2) ---
# Data preparation
df_p2_long <- df %>%
  select(PC1, Indeterminate_adj_mean, Indeterminate_adj) %>%
  rename(NullModel = Indeterminate_adj_mean, Observed = Indeterminate_adj) %>%
  pivot_longer(cols = c(NullModel, Observed), names_to = "Group", values_to = "Value") %>%
  mutate(Group = factor(Group, levels = c("NullModel", "Observed")))

# Model fitting and label generation
fit_null_p2 <- lm(Value ~ PC1, data = df_p2_long %>% filter(Group == "NullModel"))
fit_obs_p2 <- lm(Value ~ PC1, data = df_p2_long %>% filter(Group == "Observed"))

label_null_p2 <- paste0("y = ", round(coef(fit_null_p2)[1], 2), " + ", round(coef(fit_null_p2)[2], 2),
                      "x, R² = ", round(summary(fit_null_p2)$r.squared, 3), "***")
label_obs_p2 <- paste0("y = ", round(coef(fit_obs_p2)[1], 2), " + ", round(coef(fit_obs_p2)[2], 2),
                     "x, R² = ", round(summary(fit_obs_p2)$r.squared, 3), "***")

# Plotting
p2 <- ggplot(data = df_p2_long, aes(x = PC1, y = Value, color = Group, fill = Group)) +
  geom_point(alpha = 0.1, size = 1.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.3) +
  annotate("text", x = min(df$PC1, na.rm=TRUE) + 0.47 * (max(df$PC1, na.rm=TRUE) - min(df$PC1, na.rm=TRUE)), y = max(df_p2_long$Value, na.rm=TRUE), label = label_null_p2, hjust = 0, vjust = 1, size = 4, color = "black") +
  annotate("text", x = min(df$PC1, na.rm=TRUE) + 0.47 * (max(df$PC1, na.rm=TRUE) - min(df$PC1, na.rm=TRUE)), y = max(df_p2_long$Value, na.rm=TRUE) * 0.92, label = label_obs_p2, hjust = 0, vjust = 1, size = 4, color = "black") +
  scale_color_manual(values = c("NullModel" = "#707070", "Observed" = "#8B559B")) +
  scale_fill_manual(values = c("NullModel" = "#707070", "Observed" = "#8B559B")) +
  labs(x = "PC1", y = "Indeterminate_adj") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = c(0.8, 0.75),
    legend.title = element_blank(),
    legend.background = element_blank()
  ) +
  # Set Y-axis limits to match the data range
  coord_cartesian(ylim = range(df_p2_long$Value, na.rm = TRUE))

# --- Plot 3: PCA Biplot (p3) ---
# Prepare environmental data for PCA, ensuring it is numeric
env_for_pca <- df[, names(env_data_subset)]
env_for_pca <- na.omit(sapply(env_for_pca, as.numeric))

# Perform PCA
PCA_pca <- prcomp(env_for_pca, scale. = TRUE)

# Invert the direction of the PC1 axis for interpretability
PCA_pca$x[,1] <- -PCA_pca$x[,1]
PCA_pca$rotation[,1] <- -PCA_pca$rotation[,1]

# Dynamically get the percentage of variance explained by PCs
pca_summary <- summary(PCA_pca)
pc1_percent <- paste0("PC1 (", round(pca_summary$importance[2, 1] * 100, 2), "%)")
pc2_percent <- paste0("PC2 (", round(pca_summary$importance[2, 2] * 100, 2), "%)")

p3 <- fviz_pca_biplot(
  PCA_pca,
  geom = "point",
  label = "var",
  repel = TRUE,
  palette = "aaas",
  alpha.ind = 0.95,
  col.var = "black",
  col.ind = df$H_Rank,
  gradient.cols = c("#56AA5E", "#C2E6BC", "#f4f1f4", "#CAAFD4", "#8B559B"),
  ggtheme = theme_minimal(),
  pointsize = 2
) +
labs(x = pc1_percent, y = pc2_percent, color = "H_Rank") +
theme(
  legend.position = c(0.9, 0.2),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.title = element_blank(),
  plot.margin = margin(t = 0, r = 5, b = 0, l = 5),
  axis.text = element_text(size = 10),
  axis.title.y = element_text(size = 12),
  axis.title.x = element_text(size = 12),
  legend.title = element_text(size = 11),
  legend.text = element_text(size = 11)
)

# --- Combine Multiple Plots ---
# Apply common theme settings to all plots
plots <- list(p3, p2, p1)
plots <- lapply(plots, function(plot) {
  plot + theme(
    axis.title.x = element_blank(),
    plot.margin = margin(t = 0, r = 5, b = 0, l = 5),
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )
})
combined_plots <- plot_grid(plotlist = plots, ncol = 3, nrow = 1, align = 'hv')

# Create a common X-axis label
x_label <- ggdraw() +
  draw_label("PC1", x = 0.5, y = 0.5, size = 12)

# Final arrangement
final_plot <- plot_grid(
  combined_plots,
  x_label,
  ncol = 1,
  rel_heights = c(1, 0.05)
)

# Print the final plot
print(final_plot)
