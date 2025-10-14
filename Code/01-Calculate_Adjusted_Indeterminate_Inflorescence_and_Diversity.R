# =========================================================================
# Calculate Raw and Adjusted Proportions and Diversity of Inflorescences
# Last updated: 2025-10-14
# =========================================================================

# Load required libraries
library(dplyr)

# 1. Setup Working Directory and Import Data
# Set the primary working directory for the project.
setwd("D:/R_work/Inflorescence")

# Import the preprocessed dataset.
df <- read.csv(file = "Geo_Data_Trait_01_Preprocessed.csv", header = TRUE)


# 2. Calculate Initial Proportions and Diversity
# 2.1 Calculate proportions for determinate and indeterminate inflorescences
# Determinate = So_Pi + Df_Pi
df$Determinate <- df$So_Pi + df$Df_Pi

# Indeterminate = 1 - Determinate
df$Indeterminate <- 1 - df$Determinate

# 2.2 Calculate the Shannon diversity index (H)
# The formula is H = -sum(p_i * log(p_i)), where p_i is the proportion of each category.
proportion_cols <- c("So_Pi", "Df_Pi", "Co_Pi", "Sk_Pi", "Um_Pi", "Ra_Pi", "Cp_Pi", "Pa_Pi", "Ck_Pi")
df$H <- apply(df[, proportion_cols], 1, function(row) {
  # Filter out zero proportions as 0 * log(0) is defined as 0.
  p_values <- row[row > 0]
  if(length(p_values) == 0) return(0)
  -sum(p_values * log(p_values))
})


# 3. Calculate Adjusted Proportions using Empirical Bayes Shrinkage
# 3.1 Calculate the global mean proportion of indeterminate inflorescences (p0)
# This is weighted by the number of species observed in each grid cell (Species_N_Inf).
total_Indeterminate_species <- sum(df$Indeterminate * df$Species_N_Inf)
total_species <- sum(df$Species_N_Inf)
p0 <- total_Indeterminate_species / total_species

# 3.2 Calculate the mean and variance of the observed proportions
mean_p <- mean(df$Indeterminate)
var_p <- var(df$Indeterminate)

# 3.3 Estimate the alpha parameter using the method of moments
# Formula: alpha = mean_p * (1 - mean_p) / var_p - 1
alpha_optimal <- (mean_p * (1 - mean_p) / var_p - 1)
# Set a minimum threshold for alpha to prevent non-positive values.
alpha_optimal <- max(alpha_optimal, 0.01)

# 3.4 Assign the estimated alpha to be used as the weight
alpha <- alpha_optimal

# 3.5 Calculate the adjusted proportion of indeterminate inflorescences (Indeterminate_adj)
# This shrinks the observed proportions toward the global mean p0.
df$Indeterminate_adj <- (df$Indeterminate * df$Species_N_Inf + p0 * alpha) / (df$Species_N_Inf + alpha)

# 3.6 Calculate the adjusted proportion of determinate inflorescences
df$Determinate_adj <- 1 - df$Indeterminate_adj


# 4. Normalize the Shannon Diversity Index (H)
# Rank-transform H and scale it to a (0, 1) range.
df <- df %>%
  mutate(H_Rank = rank(H, ties.method = "average"),
         H_Rank_Normalized = H_Rank / (n() + 1))


# 5. Export the Processed Data
# Write the data frame with new variables to a CSV file.
write.csv(df, file = "Geo_Data_Trait_02_Adjusted.csv", row.names = FALSE)
