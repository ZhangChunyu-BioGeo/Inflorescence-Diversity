# =========================================================================
# Structural Equation Modeling (SEM) Analysis
# Last updated: 2025-10-14
# =========================================================================

# --- Load Packages ---
# lavaan is the primary package for Structural Equation Modeling in R.
library(lavaan)

# Set the primary working directory for the project.
setwd("D:/R_work/Inflorescence")

# -------------------------------------------------------------------------
# Part 1: Data Loading and Preparation
# -------------------------------------------------------------------------
# Read the comprehensive dataset containing geographic, trait, and environmental variables.
df <- read.csv("Geo_Data_Trait_and_Enveriment.csv", header = TRUE)

# Define the variables to be used in the SEM.
vars <- c("MDT", "PC1", "PD", "MCI", "Indeterminate_adj", "H_Rank")

# Normalize the selected variables to a 0-1 scale.
# This is a standard procedure to make path coefficients comparable.
normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}
df[vars] <- lapply(df[vars], normalize)

# -------------------------------------------------------------------------
# Part 2: Define and Fit the Structural Equation Model
# -------------------------------------------------------------------------
# Define the SEM model structure using lavaan syntax.
model_spec <-'
  MCI ~ PC1
  MDT ~ MCI + PC1
  PD ~ MCI + PC1
  Indeterminate_adj  ~ MDT + PC1
  H_Rank ~ PD + MCI + PC1 + Indeterminate_adj
'

# Fit the specified model to the data.
sem_fit <- sem(model_spec, data = df, estimator = "MLR")

# -------------------------------------------------------------------------
# Part 3: Display Model Results
# -------------------------------------------------------------------------
# Print a detailed summary of the model fit.
cat("--- SEM Model Summary (Standardized Coefficients) ---\n")
summary(sem_fit, standardized = TRUE, rsquare = TRUE)

# Print a comprehensive set of model fit indices.
cat("\n--- Overall Model Fit Measures ---\n")
print(fitMeasures(sem_fit))

