# ============================================================================
# Visualization of Plant-Pollinator Co-occurrence Intensity (COI/SES)
# Last updated: 2025-10-14
# ============================================================================

# --- Load Packages ---
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(circlize)

# Set the primary working directory for the project.
setwd("D:/R_work/Inflorescence")

# -------------------------------------------------------------------------
# Part 1: Load Data and Reshape
# -------------------------------------------------------------------------
# Read the Standardized Effect Size (SES/COI) results from the previous script
ses_results_df <- read.csv("Potential_Pollinators/OccNum_Plants_df_SES.csv")

# Convert the data from wide format to long format for easier plotting with ggplot2
ses_results_long <- ses_results_df %>%
  pivot_longer(
    cols = c(COI_Beetles, COI_Flies, COI_Bees, COI_Wasps, COI_Butterflies, COI_Moths),
    names_to = "Pollination",
    values_to = "co_occurrence_index",
    names_prefix = "COI_"
  )

# Set the desired order for the pollinator groups (facets in the plot)
pollinator_order <- c("Butterflies", "Bees", "Moths", "Flies", "Beetles", "Wasps")
ses_results_long <- ses_results_long %>%
   mutate(Pollination = factor(Pollination, levels = pollinator_order))

# -------------------------------------------------------------------------
# Part 2: Statistical Analysis (T-tests)
# -------------------------------------------------------------------------
# Perform t-tests to compare COI between inflorescence types for each pollinator group
ttest_results <- ses_results_long %>%
  group_by(Pollination) %>%
  do({
    # Check if there are enough data points in both groups
    if (n_distinct(.$Inflorescence) == 2 && min(table(.$Inflorescence)) > 1) {
        test <- t.test(co_occurrence_index ~ Inflorescence, data = .)
        broom::tidy(test)
    } else {
        # Return a data frame with NA if test cannot be performed
        data.frame(p.value = NA)
    }
  }) %>%
  ungroup() %>%
  mutate(
      group1 = "D",
      group2 = "I",
      p.signif = case_when(
          p.value < 0.001 ~ "***",
          p.value < 0.01  ~ "**",
          p.value < 0.05  ~ "*",
          TRUE ~ "ns"
      )
  )

# -------------------------------------------------------------------------
# Part 3: Bar Plot Visualization
# -------------------------------------------------------------------------
# Calculate summary statistics (mean, sd, n, se) for the bar plot
bar_data <- ses_results_long %>%
  group_by(Pollination, Inflorescence) %>%
  summarise(
    mean_index = mean(co_occurrence_index, na.rm = TRUE),
    sd = sd(co_occurrence_index, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  ) %>%
  mutate(se = sd / sqrt(n))

# Create the bar plot with error bars and significance annotations
p_coi_bar <- ggplot(bar_data, aes(x = Inflorescence, y = mean_index, fill = Inflorescence)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
  geom_errorbar(aes(ymin = mean_index - se, ymax = mean_index + se), 
                width = 0.2, position = position_dodge(width = 0.9)) +
  
  facet_wrap(~Pollination, nrow = 1) +
  
  scale_fill_manual(values = c("D" = "#56AA5E", "I" = "#8B559B")) + 
  
  # Add statistical significance labels
  stat_pvalue_manual(
    ttest_results, 
    label = "p.signif",
    y.position = 15, # Manually set y-position for labels
    inherit.aes = FALSE
  ) +
  
  scale_y_continuous(name = "SES of Co-occurrence Intensity", limits = c(0, 16)) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 12, face = "bold"), # Facet labels
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  ) +
  labs(fill = "Inflorescence Type")

print(p_coi_bar)


# -------------------------------------------------------------------------
# Part 4: Chord Diagram Visualization
# -------------------------------------------------------------------------
# This plot shows the strength of links between inflorescence types and pollinator groups.
cat("\nGenerating chord diagram...\n")

# Calculate the average co-occurrence index for each link
flower_traits_avg <- ses_results_long %>%
  group_by(Inflorescence, Pollination) %>%
  summarise(avg_co_occurrence = mean(co_occurrence_index, na.rm = TRUE), .groups = "drop")

# Convert the data to the matrix format required by circlize
link_matrix <- flower_traits_avg %>%
  pivot_wider(names_from = Pollination, values_from = avg_co_occurrence) %>%
  tibble::column_to_rownames(var = "Inflorescence") %>%
  as.matrix()
  
# Reorder columns to match the desired layout
link_matrix <- link_matrix[, pollinator_order]

# Set colors for the inflorescence types (grid colors)
grid_colors <- c("D" = "#56AA5E", "I" = "#8B559B")

# Open a new plotting window for the chord diagram
dev.new() 

# Reset circlize parameters before plotting
circos.clear()

# Create the chord diagram
set.seed(123) # for reproducibility of the layout
chordDiagram(
  link_matrix, 
  transparency = 0.5,
  annotationTrack = c("name", "grid"),
  grid.col = grid_colors,
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow"
)

title("Co-occurrence Intensity between Inflorescence Types and Pollinators")
cat("Chord diagram created successfully!\n")
