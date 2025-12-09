# R Script: microgastropodes_env_analyses_eduhcgalvao.R
# Author: Eduardo HC Galvao - eduardohcgalvao@gmail.com
# Version: 2.5 (22/10/2025)
# Date: 12/10/2025 (dd/mm/yyyy) #nolint
# Description: Analyses to investigate microgastropodes in three different
# pristine river basins in Brazil.

# The workflow includes:
# 1) Data import and cleaning
# 2) Environmental data analyses: test differences between river basins
# (pH, temperature, salinity, dissolved oxygen, suspended solids,
# sediment particle size) in those three river basins and across seasons.
# 3) Diversity analyses: compare differences in shannon-wiener, number of
# species, margalef richness and pielou evenness between river basins
# and seasons.
# 4) Community composition: test if the community composition differs 
# between river basins and seasons using PERMANOVA and NMDS.
# 5) Relationship between environmental factors with shannon-wiener and number
# of species using GAMs.
# 6) Relationship between environmental factors and community composition
# using distance-based redundancy analysis (dbRDA).
# 7) Investigate which species contribute the most for the differences between
# river basins and seasons using SIMPER.

# Input data:
# 1) Environmental factors data: "env_data_microgastropodes.csv"
# 2) Community data: "community_data_microgastropodes.csv"
# 3) Macronutrients data: "macronutrients_data_microgastropodes.csv"

# If one the following packages are not installed run:
# install.packages("package") #nolint
library(ggplot2) # For plotting
library(dplyr) # For data manipulation
library(tidyr) # For data reshaping
library(reshape2) # For data manipulation
library(car) # For VIF analysis and additional statistical tests
library(vegan)     # For dbRDA, PERMANOVA, diversity indices
library(cowplot)   # For combining plots
library(mgcv) # For GAM models

set.seed(123) # For reproducibility

###############################################################################
############################ Load and prepare data ############################
###############################################################################

# Set working directory
setwd("C:/Users/Kelmo/Desktop/eduardo/Microgastropodes_as_environmental_indicators")

# Load data
env_var_data <- read.csv("data/env_data_microgastropodes.csv")
community_data <- read.csv("data/community_data_microgastropodes.csv")

# Get species names and abundance data
species_names <- community_data[,1]
abundance_data <- community_data[,-1]

# Create sampling point information
sampling_points <- colnames(abundance_data)
site_codes <- c("C5", "C6", "C7", "C8", "C9", "C10", "N1", "N2",
                        "N3", "N4", "S1", "S2", "S3", "S4", "S5")
season_codes <- substr(sampling_points, nchar(sampling_points),
                                            nchar(sampling_points))
site_names <- substr(sampling_points, 1, nchar(sampling_points)-1)
seasons <- ifelse(season_codes == "A", "Wet", "Dry")

# Create a community matrix for analyses
community_matrix <- t(abundance_data)
colnames(community_matrix) <- species_names
rownames(community_matrix) <- sampling_points

# Create site-season dataframe
site_season_data <- data.frame(
  sampling_points = sampling_points,
  site = site_names,
  season = seasons,
  stringsAsFactors = FALSE
)

# Add environment mapping for macronutrients merge
# Map sites to environments based on the data structure
site_to_environment <- data.frame(
  site = c("N1", "N2", "N3", "N4", "C5", "C6", "C7", "C8",
                    "C9", "C10", "S1", "S2", "S3", "S4", "S5"),
  environment = c(rep("Serinhaem", 4), rep("Sorojo", 6), rep("Marau", 5)),
  stringsAsFactors = FALSE
)

# Add environment to site_season_data
site_season_data <- merge(site_season_data,
                          site_to_environment,
                          by = "site",
                          all.x = TRUE)

# Merge environmental variables data with site-season information
combined_data <- merge(site_season_data,
                       env_var_data, 
                       by.x = c("site", "season"),
                       by.y = c("sampling_points", "season"), 
                       suffixes = c("", "_env_var"))

# Ensure the order matches the community matrix
combined_data <- combined_data[match(site_season_data$sampling_points,
                                        combined_data$sampling_points),]

# Create environmental matrix for analyses
env_matrix <- combined_data[, c("temperature",
                                "salinity",
                                "pH",
                                "dissolved_o2",
                                "suspended_solids", 
                                "sediment_particle_size")]
rownames(env_matrix) <- combined_data$sampling_points

# Write down a summary of abundance data for season and river basin
abundance_summary <- data.frame(
  sampling_points = combined_data$sampling_points,
  environment = combined_data$environment,
  season = combined_data$season,
  total_abundance = rowSums(community_matrix[combined_data$sampling_points, ])
)

###############################################################################
######################### Environmental data analyses #########################
###############################################################################

# Calculate mean and sd of environmental variables by environment
environments <- c("Serinhaem", "Sorojo", "Marau")
seasons <- c("Dry", "Wet")
env_variables <- c("temperature", "salinity", "pH", "dissolved_o2",
                   "suspended_solids", "sediment_particle_size")

# Define function to compute statistics for groups
compute_stats <- function(data, group_var, group_values, var_list) {
  result <- data.frame()
  
  for (val in group_values) {
    subset_data <- data[data[[group_var]] == val, ]
    
    for (var in var_list) {
      if (var %in% colnames(subset_data)) {
        var_values <- subset_data[[var]]
        var_values <- var_values[!is.na(var_values)]
        
        if (length(var_values) > 0) {
          result <- rbind(result, data.frame(
            Group = val,
            Variable = var,
            Mean = round(mean(var_values), 3),
            SD = round(sd(var_values), 3),
            N = length(var_values)
          ))
        }
      }
    }
  }
  
  return(result)
}

################################################################################

# Remove P column because of constant values
combined_data_filtered <- combined_data[, !names(combined_data) %in% "P"]

# Define function to perform statistical tests
perform_statistical_tests <- function(data, variables, grouping_var) {
  result <- data.frame(
    Variable = character(),
    Shapiro_P = character(),
    Test = character(),
    Statistic = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (var in variables) {
    if (var %in% colnames(data)) {
      var_values <- data[[var]]
      group_values <- data[[grouping_var]]

      
      # Check if we have sufficient data and variation
      if (length(var_values) > 2 && length(unique(group_values)) > 1) {
        tryCatch({
          # First fit ANOVA to get residuals for normality testing
          anova_model <- aov(var_values ~ group_values)
          anova_residuals <- residuals(anova_model)
          
          # Test normality of residuals
          if (length(anova_residuals) >= 3 && length(anova_residuals) <= 5000) {
            shapiro_p <- shapiro.test(anova_residuals)$p.value
            is_normal <- shapiro_p > 0.05
          } else {
            shapiro_p <- NA  # Shapiro test not applicable
            is_normal <- FALSE
          }
          
          # Choose appropriate test based on residual normality
          if (is_normal) {
            # Use ANOVA results (already computed above)
            anova_summary <- summary(anova_model)
            p_value <- anova_summary[[1]][["Pr(>F)"]][1]
            statistic <- anova_summary[[1]][["F value"]][1]
            test_type <- "ANOVA"
          } else {
            # Use Kruskal-Wallis for non-normal residuals
            kruskal_test <- kruskal.test(var_values ~ group_values)
            p_value <- kruskal_test$p.value
            statistic <- kruskal_test$statistic
            test_type <- "Kruskal-Wallis"
          }
          
          result <- rbind(result, data.frame(
            Variable = var,
            Shapiro_P = ifelse(is.na(shapiro_p), "NA", 
                              ifelse(shapiro_p < 0.001, 
                                    format(shapiro_p,
                                           scientific = TRUE,
                                           digits = 3),
                                    as.character(round(shapiro_p, 4)))),
            Test = test_type,
            Statistic = round(statistic, 4),
            P_value = ifelse(p_value < 0.001,
                            format(p_value, scientific = TRUE, digits = 3),
                            round(p_value, 4))
          ))
        }, error = function(e) {
          cat("Error testing", var, ":", e$message, "\n")
        })
      }
    }
  }
  
  return(result)
}

# Test differences in environmental variables between environments
env_dif_results <- perform_statistical_tests(combined_data,
                                             env_variables,
                                             "environment")

# Test differences in environmental variables between seasons
env_season_dif_results <- perform_statistical_tests(combined_data,
                                                    env_variables,
                                                    "season")

#print(env_dif_results)
#print(env_season_dif_results)

# Save results
#write.csv(env_dif_results, "results/environmental_variables_river_anova_kruskal_results.csv", row.names = FALSE)
#write.csv(env_season_dif_results, "results/environmental_variables_season_anova_kruskal_results.csv", row.names = FALSE)

# Two way ANOVA (season and environment) for the environmental variables
two_way_anova_results <- data.frame(
  Variable = character(),
  Season_P = numeric(),
  Environment_P = numeric(),
  Interaction_P = numeric(),
  stringsAsFactors = FALSE
)

for (var in env_variables) {
  if (var %in% colnames(combined_data)) {
    var_values <- combined_data[[var]]
    season_values <- combined_data[["season"]]
    environment_values <- combined_data[["environment"]]
  
    # Check if we have sufficient data and variation
    if (length(var_values) > 2 && length(unique(season_values)) > 1 && length(unique(environment_values)) > 1) {
      tryCatch({
        # Fit two-way ANOVA model
        anova_model <- aov(var_values ~ season_values * environment_values)
        anova_summary <- summary(anova_model)
        
        # Extract p-values
        season_p <- anova_summary[[1]][["Pr(>F)"]][1]
        environment_p <- anova_summary[[1]][["Pr(>F)"]][2]
        interaction_p <- anova_summary[[1]][["Pr(>F)"]][3]
        
        two_way_anova_results <- rbind(two_way_anova_results, data.frame(
          Variable = var,
          Season_P = ifelse(season_p < 0.001,
                            format(season_p, scientific = TRUE, digits = 3),
                            round(season_p, 4)),
          Environment_P = ifelse(environment_p < 0.001,
                                 format(environment_p, scientific = TRUE, digits = 3),
                                 round(environment_p, 4)),
          Interaction_P = ifelse(interaction_p < 0.001,
                                 format(interaction_p, scientific = TRUE, digits = 3),
                                 round(interaction_p, 4))
        ))
      }, error = function(e) {
        cat("Error in two-way ANOVA for", var, ":", e$message, "\n")
      })
    }
  }
}

#print(two_way_anova_results)

# Save results
#write.csv(two_way_anova_results, "results/environmental_variables_two_way_anova_results.csv", row.names = FALSE)

# Create boxplot panel for environmental variables by environment and season

# Function to calculate boxplot statistics
calculate_boxplot_stats <- function(data, group_col, value_col) {
  data %>%
    group_by(!!sym(group_col)) %>%
    summarise(
      ymin = min(!!sym(value_col), na.rm = TRUE),
      lower = quantile(!!sym(value_col), 0.25, na.rm = TRUE),
      middle = median(!!sym(value_col), na.rm = TRUE),
      upper = quantile(!!sym(value_col), 0.75, na.rm = TRUE),
      ymax = max(!!sym(value_col), na.rm = TRUE),
      .groups = "drop"
    )
}

# Function to format p-values
format_p_value <- function(p_value) {
  if (is.numeric(p_value) && !is.na(p_value)) {
    if (p_value < 0.001) {
      return("p < 0.001")
    } else {
      return(paste("p =", round(p_value, 3)))
    }
  } else {
    return("p = N/A")
  }
}

# Create custom Y-axis labels for environmental variables
env_var_labels <- c(
  "temperature" = "Temperature (°C)",
  "salinity" = "Salinity",
  "pH" = "pH", 
  "dissolved_o2" = "Dissolved O₂ (mg/L)",
  "suspended_solids" = "Suspended solids (mg/L)",
  "sediment_particle_size" = "Sediment particle size (φ)"
)

# Panel labels A-F
panel_labels <- c("A", "B", "C", "D", "E", "F")

# Prepare data for plotting
plot_data <- combined_data %>%
  select(environment, season, all_of(env_variables)) %>%
  pivot_longer(cols = all_of(env_variables), 
               names_to = "Variable", 
               values_to = "Value")

# Create individual plots for each environmental variable
plot_list <- list()

for (i in 1:length(env_variables)) {
  var_name <- env_variables[i]
  
  # Filter data for current variable
  var_data <- plot_data %>% filter(Variable == var_name)
  
  # Create interaction variable for grouping
  var_data$interaction <- paste(var_data$environment, var_data$season, sep = "_")
  
  # Calculate boxplot statistics
  boxplot_stats <- calculate_boxplot_stats(var_data, "interaction", "Value")
  boxplot_stats$environment <- sub("_.*", "", boxplot_stats$interaction)
  boxplot_stats$season <- sub(".*_", "", boxplot_stats$interaction)
  
  # Get interaction p-value from two-way ANOVA results
  interaction_p_text <- ""
  
  if (exists("two_way_anova_results") && nrow(two_way_anova_results) > 0) {
    var_results <- two_way_anova_results[two_way_anova_results$Variable == var_name, ]
    if (nrow(var_results) > 0) {
      # Extract interaction p-value and convert to numeric if needed
      interaction_p <- var_results$Interaction_P
      
      # Convert to numeric if in character format
      if (is.character(interaction_p)) {
        interaction_p <- as.numeric(gsub("[^0-9.e-]", "", interaction_p))
      }
      
      # Format interaction p-value
      interaction_p_text <- paste(format_p_value(interaction_p))
    }
  }
  
  # Create plot
  p <- ggplot(var_data, aes(x = environment, fill = season)) +
    geom_boxplot(
      data = boxplot_stats,
      aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
      stat = "identity", alpha = 0.8, outlier.shape = 16, outlier.size = 2,
      position = position_dodge(width = 0.8)
    ) +
    geom_jitter(aes(y = Value), alpha = 0.6, size = 2,
                position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +
    scale_fill_manual(values = c("Wet" = "#922ecc", "Dry" = "#3873d1"),
                      name = "Season") +
    labs(x = "River Basin",
         y = env_var_labels[var_name]) +
    annotate("text", x = Inf, y = Inf, label = interaction_p_text, face = "bold",
             hjust = 1.05, vjust = 2.5, size = 5, color = "black") +
    annotate("text", x = -Inf, y = Inf, label = panel_labels[i], face = "bold",
             hjust = -0.5, vjust = 1.5, size = 6, color = "black") +
    theme_classic() +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 14, face = "bold"),
      axis.title.x = if(i > 3) element_text() else element_blank(),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  plot_list[[i]] <- p
}

# Arrange plots using cowplot
env_boxplot_panel <- plot_grid(
  plot_list[[1]], plot_list[[2]], plot_list[[3]],
  plot_list[[4]], plot_list[[5]], plot_list[[6]],
  nrow = 2, ncol = 3,
  align = "hv",
  axis = "tb"
)

# Extract legend from one of the plots
legend <- get_legend(
  plot_list[[1]] + 
    theme(legend.position = "bottom",
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12))
)

# Combine plots with legend at bottom
env_boxplot_panel <- plot_grid(
  env_boxplot_panel,
  legend,
  nrow = 2,
  rel_heights = c(1, 0.1)
)

# Save the panel with interaction p-values
#ggsave("results/environmental_variables_boxplot_panel.png", env_boxplot_panel, width = 15, height = 10, dpi = 1200)

##############################################################################
######################### Diversity Indices Analysis #########################
##############################################################################

# Calculate diversity indices
# Filter to complete environmental cases first
env_complete <- na.omit(env_matrix)
community_complete <- community_matrix[rownames(env_complete), ]

# Calculate diversity indices on the filtered data
shannon_wiener <- diversity(community_complete, index = "shannon")
number_of_species <- specnumber(community_complete)
margalef_richness <- (number_of_species - 1) / log(rowSums(community_complete))
pielou_evenness <- shannon_wiener / log(number_of_species)

# Get season and river basin vectors for filtered samples
season_vector <- combined_data$season[match(rownames(community_complete), combined_data$sampling_points)]
river_basin_vector <- combined_data$environment[match(rownames(community_complete), combined_data$sampling_points)]

diversity_data <- data.frame(
  sampling_points = rownames(community_complete),
  season = season_vector,
  river_basin = river_basin_vector,
  shannon = shannon_wiener,
  number_of_species = number_of_species,
  evenness = pielou_evenness,
  margalef = margalef_richness
)

# Calculate mean and standard deviation for diversity indices by season
diversity_summary_season <- diversity_data %>%
  group_by(season) %>%
  summarise(
    mean_shannon = mean(shannon, na.rm = TRUE),
    sd_shannon = sd(shannon, na.rm = TRUE),
    mean_number_of_species = mean(number_of_species, na.rm = TRUE),
    sd_number_of_species = sd(number_of_species, na.rm = TRUE),
    mean_evenness = mean(evenness, na.rm = TRUE),
    sd_evenness = sd(evenness, na.rm = TRUE),
    mean_margalef = mean(margalef, na.rm = TRUE),
    sd_margalef = sd(margalef, na.rm = TRUE),
    n = n()
  )

# Calculate mean and standard deviation for diversity indices by river basin
diversity_summary_river_basin <- diversity_data %>%
  group_by(river_basin) %>%
  summarise(
    mean_shannon = mean(shannon, na.rm = TRUE),
    sd_shannon = sd(shannon, na.rm = TRUE),
    mean_number_of_species = mean(number_of_species, na.rm = TRUE),
    sd_number_of_species = sd(number_of_species, na.rm = TRUE),
    mean_evenness = mean(evenness, na.rm = TRUE),
    sd_evenness = sd(evenness, na.rm = TRUE),
    mean_margalef = mean(margalef, na.rm = TRUE),
    sd_margalef = sd(margalef, na.rm = TRUE),
    n = n()
  )

# Calculate mean and standard deviation for diversity indices by season × river basin
diversity_summary_season_river_basin <- diversity_data %>%
  group_by(season, river_basin) %>%
  summarise(
    mean_shannon = mean(shannon, na.rm = TRUE),
    sd_shannon = sd(shannon, na.rm = TRUE),
    mean_number_of_species = mean(number_of_species, na.rm = TRUE),
    sd_number_of_species = sd(number_of_species, na.rm = TRUE),
    mean_evenness = mean(evenness, na.rm = TRUE),
    sd_evenness = sd(evenness, na.rm = TRUE),
    mean_margalef = mean(margalef, na.rm = TRUE),
    sd_margalef = sd(margalef, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

# Save all diversity summaries
#write.csv(diversity_data, "results/diversity_indices_data.csv", row.names = FALSE)
#write.csv(diversity_summary_season, "results/diversity_summary_by_season.csv", row.names = FALSE)
#write.csv(diversity_summary_river_basin, "results/diversity_summary_by_river_basin.csv", row.names = FALSE)
#write.csv(diversity_summary_season_river_basin, "results/diversity_summary_by_season_river_basin.csv", row.names = FALSE)

# Test differences in diversity indices between seasons and river basins
diversity_variables <- c("shannon", "number_of_species", "evenness", "margalef")
diversity_formatted_season <- perform_statistical_tests(diversity_data,
                                                         diversity_variables,
                                                         "season")
diversity_formatted_river_basin <- perform_statistical_tests(diversity_data,
                                                            diversity_variables,
                                                            "river_basin")

# Two way ANOVA (season and river basin) for diversity indices
diversity_two_way_anova_results <- data.frame(
  Variable = character(),
  Season_P = numeric(),
  River_Basin_P = numeric(),
  Interaction_P = numeric(),
  stringsAsFactors = FALSE
)

for (var in diversity_variables) {
  if (var %in% colnames(diversity_data)) {
    var_values <- diversity_data[[var]]
    season_values <- diversity_data[["season"]]
    river_basin_values <- diversity_data[["river_basin"]]
    
    # Check if we have sufficient data and variation
    if (length(var_values) > 2 && length(unique(season_values)) > 1 && length(unique(river_basin_values)) > 1) {
      tryCatch({
        # Fit two-way ANOVA model
        anova_model <- aov(var_values ~ season_values * river_basin_values)
        anova_summary <- summary(anova_model)
        
        # Extract p-values
        season_p <- anova_summary[[1]][["Pr(>F)"]][1]
        river_basin_p <- anova_summary[[1]][["Pr(>F)"]][2]
        interaction_p <- anova_summary[[1]][["Pr(>F)"]][3]
        
        diversity_two_way_anova_results <- rbind(diversity_two_way_anova_results, data.frame(
          Variable = var,
          Season_P = ifelse(season_p < 0.001,
                            format(season_p, scientific = TRUE, digits = 3),
                            round(season_p, 4)),
          River_Basin_P = ifelse(river_basin_p < 0.001,
                                 format(river_basin_p, scientific = TRUE, digits = 3),
                                 round(river_basin_p, 4)),
          Interaction_P = ifelse(interaction_p < 0.001,
                                 format(interaction_p, scientific = TRUE, digits = 3),
                                 round(interaction_p, 4))
        ))
      }, error = function(e) {
        cat("Error in two-way ANOVA for", var, ":", e$message, "\n")
      })
    }
  }
}

#print(diversity_two_way_anova_results)

# Save results
#write.csv(diversity_formatted_season, "results/diversity_indices_season_anova_kruskal_results.csv", row.names = FALSE)
#write.csv(diversity_formatted_river_basin, "results/diversity_indices_river_basin_anova_kruskal_results.csv", row.names = FALSE)
#write.csv(diversity_two_way_anova_results, "results/diversity_indices_two_way_anova_results.csv", row.names = FALSE)

# Create boxplot panel for diversity indices by river basin and season

# Create custom Y-axis labels for diversity indices
diversity_var_labels <- c(
  "shannon" = "Shannon Diversity",
  "number_of_species" = "Number of Species",
  "evenness" = "Pielou's Evenness",
  "margalef" = "Margalef Richness"
)

# Panel labels A-D
diversity_panel_labels <- c("A", "B", "C", "D")

# Prepare data for plotting
diversity_plot_data <- diversity_data %>%
  select(river_basin, season, all_of(diversity_variables)) %>%
  pivot_longer(cols = all_of(diversity_variables), 
               names_to = "Variable", 
               values_to = "Value")

# Create individual plots for each diversity index
diversity_plot_list <- list()

for (i in 1:length(diversity_variables)) {
  var_name <- diversity_variables[i]
  
  # Filter data for current variable
  var_data <- diversity_plot_data %>% filter(Variable == var_name)
  
  # Create interaction variable for grouping
  var_data$interaction <- paste(var_data$river_basin, var_data$season, sep = "_")
  
  # Calculate boxplot statistics
  boxplot_stats <- calculate_boxplot_stats(var_data, "interaction", "Value")
  boxplot_stats$river_basin <- sub("_.*", "", boxplot_stats$interaction)
  boxplot_stats$season <- sub(".*_", "", boxplot_stats$interaction)
  
  # Get interaction p-value from two-way ANOVA results
  p_value_text <- ""
  
  if (exists("diversity_two_way_anova_results") && nrow(diversity_two_way_anova_results) > 0) {
    var_results <- diversity_two_way_anova_results[diversity_two_way_anova_results$Variable == var_name, ]
    if (nrow(var_results) > 0) {
      # Extract interaction p-value and convert to numeric if needed
      interaction_p <- var_results$Interaction_P
      
      # Convert to numeric if in character format
      if (is.character(interaction_p)) {
        interaction_p <- as.numeric(gsub("[^0-9.e-]", "", interaction_p))
      }
      
      # Format interaction p-value (without "Interaction:" prefix)
      p_value_text <- format_p_value(interaction_p)
    } else {
      p_value_text <- "p = NS"
    }
  } else {
    p_value_text <- "p = NS"
  }
  
  # Create plot
  p <- ggplot(var_data, aes(x = river_basin, fill = season)) +
    geom_boxplot(
      data = boxplot_stats,
      aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
      stat = "identity", alpha = 0.8, outlier.shape = 16, outlier.size = 2,
      position = position_dodge(width = 0.8)
    ) +
    geom_jitter(aes(y = Value), alpha = 0.6, size = 2,
                position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +
    scale_fill_manual(values = c("Wet" = "#922ecc", "Dry" = "#3873d1"),
                      name = "Season") +
    labs(x = "River Basin",
         y = diversity_var_labels[var_name]) +
    annotate("text", x = Inf, y = Inf, label = p_value_text, face = "bold",
             hjust = 1.05, vjust = 2.5, size = 5, color = "black") +
    annotate("text", x = -Inf, y = Inf, label = diversity_panel_labels[i], face = "bold",
             hjust = -0.5, vjust = 1.5, size = 6, color = "black") +
    theme_classic() +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 14, face = "bold"),
      axis.title.x = if(i > 2) element_text() else element_blank(),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  diversity_plot_list[[i]] <- p
}

# Arrange plots using cowplot
diversity_boxplot_panel <- plot_grid(
  diversity_plot_list[[1]], diversity_plot_list[[2]], 
  diversity_plot_list[[3]], diversity_plot_list[[4]],
  nrow = 2, ncol = 2,
  align = "hv",
  axis = "tb"
)

# Extract legend from one of the plots
diversity_legend <- get_legend(
  diversity_plot_list[[1]] + 
    theme(legend.position = "bottom",
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12))
)

# Combine plots with legend at bottom
diversity_boxplot_panel <- plot_grid(
  diversity_boxplot_panel,
  diversity_legend,
  nrow = 2,
  rel_heights = c(1, 0.1)
)

# Save the panel
#ggsave("results/diversity_indices_boxplot_panel.png", diversity_boxplot_panel, width = 12, height = 10, dpi = 1200)

################################################################################
######################### Community Composition Analysis #######################
################################################################################

# Transform community data to fourth-root
env_complete <- na.omit(env_matrix)
community_complete <- community_matrix[rownames(env_complete), ]
community_transformed <- community_complete^0.25

# Prepare data for PERMANOVA
# Create vectors that match the filtered community data
season_vector <- combined_data$season[match(rownames(community_complete),
                                              combined_data$sampling_points)]
site_vector <- combined_data$site[match(rownames(community_complete),
                                              combined_data$sampling_points)]
environment_vector <- combined_data$environment[match(rownames(community_complete),
                                                        combined_data$sampling_points)]

# Create distance matrix
community_dist <- vegdist(community_transformed, method = "bray")

# PERMANOVA testing season effects
permanova_season <- adonis2(community_dist ~ season_vector, permutations = 9999)
#print(permanova_season)

# PERMANOVA testing environment type effects
permanova_environment <- adonis2(community_dist ~ environment_vector,
                                                          permutations = 9999)
#print(permanova_environment)

# Combined PERMANOVA with season and environment
permanova_combined <- adonis2(community_dist ~ season_vector + environment_vector,
                                                                permutations = 9999)
#print(permanova_combined)

# Test interaction between season and environment
permanova_interaction <- adonis2(community_dist ~ season_vector * environment_vector,
                                                                   permutations = 9999)
#print(permanova_interaction)

# PERMADISPER to test homogeneity of dispersions
betadisper_season <- betadisper(community_dist, season_vector)
permutest_season <- permutest(betadisper_season, permutations = 9999)
#print(permutest_season)

# PERMADISPER to test homogeneity of dispersions
betadisper_environment <- betadisper(community_dist, environment_vector)
permutest_environment <- permutest(betadisper_environment, permutations = 9999)
#print(permutest_environment)

# Save PERMANOVA results to CSV
#write.csv(as.data.frame(permanova_season), "results/permanova_season_results.csv", row.names = TRUE)
#write.csv(as.data.frame(permanova_environment), "results/permanova_environmental_results.csv", row.names = TRUE)
#write.csv(as.data.frame(permanova_combined), "results/permanova_combined_results.csv", row.names = TRUE)
#write.csv(as.data.frame(permanova_interaction), "results/permanova_interaction_results.csv", row.names = TRUE)

# Save PERMADISPER results to CSV
#write.csv(as.data.frame(permutest_season$tab), "results/permadisper_season_results.csv", row.names = TRUE)
#write.csv(as.data.frame(permutest_environment$tab), "results/permadisper_environment_results.csv", row.names = TRUE)

################################################################################
############################### NMDS Plotting ##################################
################################################################################

# Perform NMDS with Bray-Curtis dissimilarity
nmds_result <- metaMDS(community_transformed, distance = "bray", k = 2, trymax = 100)
#print(nmds_result)

# Create NMDS plot with environmental vectors
# Extract NMDS scores for plotting
nmds_site_scores <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_species_scores <- as.data.frame(scores(nmds_result, display = "species"))

# Add labels and season information
nmds_site_scores$site_label <- rownames(nmds_site_scores)
nmds_site_scores$season <- season_vector
nmds_site_scores$site <- site_vector
nmds_site_scores$environment <- environment_vector

# Create NMDS plot with environmental vectors and PERMANOVA results
# Extract PERMANOVA combined results for subtitle
permanova_r2 <- round(permanova_combined$R2[1], 3)  # Total R² explained
permanova_p_raw <- permanova_combined$`Pr(>F)`[1]       # Raw p-value without rounding

# Format p-value for display (use raw p-value for comparison)
permanova_p_text <- if(permanova_p_raw < 0.001) "p < 0.001" else paste("p =", round(permanova_p_raw, 3))

nmds_plot <- ggplot() +
  geom_point(data = nmds_site_scores, aes(x = NMDS1, y = NMDS2, 
                                          color = environment, shape = season),
             size = 4, alpha = 0.8, stroke = 1.5) +
  stat_ellipse(data = nmds_site_scores, aes(x = NMDS1, y = NMDS2, color = environment),
               level = 0.95, linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_text(data = nmds_site_scores, aes(x = NMDS1, y = NMDS2, label = site_label),
            size = 3, vjust = -0.8, hjust = 0.5, color = "black") +
  scale_color_manual(values = c("Serinhaem" = "#18b798", "Sorojo" = "#1296ff", "Marau" = "#da4ef3"),
                     name = "River Basin") +
  scale_shape_manual(values = c("Wet" = 16, "Dry" = 17),
                     name = "Season") +
  labs(title = "NMDS Ordination",
       subtitle = paste0("Stress = ", round(nmds_result$stress, 3), 
                        " | PERMANOVA: R² = ", permanova_r2, ", ", permanova_p_text),
       x = "NMDS1", y = "NMDS2") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    legend.box = "horizontal"
  ) +
  # Add origin lines
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "gray50") +
  guides(color = guide_legend(override.aes = list(size = 5)),
         shape = guide_legend(override.aes = list(size = 5)))

# Save the plot
#ggsave("results/nmds_plot_ggplot_PERMANOVA.png", nmds_plot, width = 12, height = 8, dpi = 1200)

################################################################################
################################# VIF Analysis #################################
################################################################################

# Check multicollinearity of environmental variables using VIF
# VIF > 7 indicates multicollinearity

# Calculate VIF for each environmental variable
vif_results <- data.frame(Variable = character(), VIF = numeric(), stringsAsFactors = FALSE)

for (i in 1:ncol(env_complete)) {
  # Create a model with all other variables as predictors
  other_vars <- env_complete[, -i, drop = FALSE]
  if (ncol(other_vars) > 0 && nrow(other_vars) > ncol(other_vars)) {
    target_var <- env_complete[, i]
    vif_model <- lm(target_var ~ ., data = as.data.frame(other_vars))
    r_squared <- summary(vif_model)$r.squared
    vif_value <- 1 / (1 - r_squared)
    vif_results <- rbind(vif_results, data.frame(Variable = colnames(env_complete)[i], VIF = vif_value))
  }
}

#print(vif_results)

# Save VIF results to CSV
#write.csv(vif_results, "results/vif_results.csv", row.names = FALSE)

# Filter variables with VIF < 7
low_vif_vars <- vif_results$Variable[vif_results$VIF < 7]
env_filtered <- env_complete[, low_vif_vars, drop = FALSE]

# Also remove P from the env_filtered
env_filtered <- env_filtered[, !colnames(env_filtered) %in% "P", drop = FALSE]

################################################################################
############################# Diversity Modeling ###############################
################################################################################

# Prepare model data
model_data <- data.frame(
  shannon = diversity_data$shannon,
  number_of_species = diversity_data$number_of_species,
  temperature = env_filtered$temperature,
  salinity = env_filtered$salinity,
  dissolved_o2 = env_filtered$dissolved_o2,
  suspended_solids = env_filtered$suspended_solids,
  sediment_particle_size = env_filtered$sediment_particle_size,
  pH = env_complete$pH
)

# Function to fit GAM and create plot for one variable and one response
make_gam_plot <- function(response, predictor, x_label, y_label) {
  # Build formula dynamically
  formula <- as.formula(paste(response, "~ s(", predictor, ")", sep = ""))
  
  # Fit GAM model
  model <- gam(formula, data = model_data)
  stats <- summary(model)
  
  # Extract p-value and explained deviance
  p_value <- round(stats$s.table[1, "p-value"], 3)
  deviance <- round(stats$dev.expl * 100, 1)
  p_text <- ifelse(p_value < 0.001, "p < 0.001", paste("p =", p_value))
  annotation <- paste0(p_text, "\nDev. = ", deviance, "%")
  
  # Create plot
  plot <- ggplot(model_data, aes_string(x = predictor, y = response)) +
    geom_point(size = 3) +
    geom_smooth(method = "gam", colour = "red", linewidth = 3) +
    labs(x = x_label, y = y_label) +
    annotate("text", x = Inf, y = Inf, label = annotation, 
             hjust = 1.1, vjust = 1.1, size = 5, fontface = "bold") +
    theme_classic() +
    theme(
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 16, face = "bold"),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  
  return(plot)
}

# Function to generate combined Shannon/Number of Species plots for each variable
make_variable_panel <- function(var, x_label) {
  shannon_plot <- make_gam_plot("shannon", var, x_label, "Shannon Diversity")
  richness_plot <- make_gam_plot("number_of_species", var, x_label, "Number of Species")

  combined <- plot_grid(
    shannon_plot + theme(axis.title.x = element_blank()),
    richness_plot,
    ncol = 1, align = "hv", label_size = 16, label_fontface = "bold"
  )
  return(combined)
}

# Generate all variable panels
plots <- list(
  DissolvedO2 = make_variable_panel("dissolved_o2", "Dissolved O₂ (mg/L)"),
  Temperature = make_variable_panel("temperature", "Temperature (°C)"),
  Salinity = make_variable_panel("salinity", "Salinity (ppt)"),
  SuspSolids = make_variable_panel("suspended_solids", "Suspended Solids (mg/L)"),
  Sediment = make_variable_panel("sediment_particle_size", "Sediment Particle Size (φ)"),
  pH = make_variable_panel("pH", "pH")
)

# Combine all panels into one final figure
all_gam_combined_plot <- plot_grid(
  plots$DissolvedO2, plots$Temperature, plots$Salinity,
  plots$SuspSolids, plots$Sediment, plots$pH,
  ncol = 6, align = "hv", label_size = 16, label_fontface = "bold"
)

# Save figure
#ggsave("results/all_gam_models.png", all_gam_combined_plot, width = 18, height = 10, dpi = 1200)

# Function to extract GAM model statistics
extract_gam_stats <- function(response, predictor, predictor_name) {
  # Build formula dynamically
  formula <- as.formula(paste(response, "~ s(", predictor, ")", sep = ""))
  
  # Fit GAM model
  model <- gam(formula, data = model_data)
  stats <- summary(model)
  
  # Extract statistics
  p_value <- stats$s.table[1, "p-value"]
  deviance <- stats$dev.expl * 100
  edf <- stats$s.table[1, "edf"]  # Effective degrees of freedom
  f_value <- stats$s.table[1, "F"]
  aic <- AIC(model)
  
  # Create results data frame
  results <- data.frame(
    Response = response,
    Predictor = predictor_name,
    P_value = p_value,
    Deviance_Explained_Percent = round(deviance, 2),
    Effective_DF = round(edf, 3),
    F_value = round(f_value, 3),
    AIC = round(aic, 2),
    Significance = ifelse(p_value < 0.001, "***", 
                         ifelse(p_value < 0.01, "**", 
                               ifelse(p_value < 0.05, "*", "ns")))
  )
  
  return(results)
}

# Extract GAM results for all variable combinations
gam_results <- data.frame()

# List of predictors and their names
predictors <- list(
  list(var = "dissolved_o2", name = "Dissolved O2"),
  list(var = "temperature", name = "Temperature"),
  list(var = "salinity", name = "Salinity"),
  list(var = "suspended_solids", name = "Suspended Solids"),
  list(var = "sediment_particle_size", name = "Sediment Particle Size"),
  list(var = "pH", name = "pH")
)

# Response variables
responses <- c("shannon", "number_of_species")

# Extract statistics for all combinations
for (pred in predictors) {
  for (resp in responses) {
    tryCatch({
      result <- extract_gam_stats(resp, pred$var, pred$name)
      gam_results <- rbind(gam_results, result)
    }, error = function(e) {
      cat("Error fitting GAM for", resp, "~", pred$name, ":", e$message, "\n")
    })
  }
}

# Save GAM results to CSV
#write.csv(gam_results, "results/gam_model_results.csv", row.names = FALSE)

################################################################################
################## Distance-based Redundancy Analysis (dbRDA) ##################
################################################################################

# Perform dbRDA with VIF-filtered environmental variables
dbrda_result <- dbrda(community_dist ~ ., data = as.data.frame(env_filtered))
#print(summary(dbrda_result))

# Test significance of dbRDA axes
dbrda_significance <- anova(dbrda_result, permutations = 9999)
#print(dbrda_significance)

# Test significance of individual environmental variables
dbrda_var_significance <- anova(dbrda_result, by = "terms", permutations = 9999)
#print(dbrda_var_significance)

# Save dbRDA results to CSV
dbrda_summary <- summary(dbrda_result)
#write.csv(as.data.frame(dbrda_summary$cont$importance), "results/dbrda_eigenvalues.csv", row.names = TRUE)
#write.csv(as.data.frame(dbrda_var_significance), "results/dbrda_variable_significance.csv", row.names = TRUE)
#write.csv(as.data.frame(dbrda_significance), "results/dbrda_overall_significance.csv", row.names = TRUE)

# Create dbRDA biplot
# Extract dbRDA scores for plotting
site_scores <- as.data.frame(scores(dbrda_result, display = "sites"))
species_scores <- as.data.frame(scores(dbrda_result, display = "species"))
env_scores <- as.data.frame(scores(dbrda_result, display = "bp"))

# Add season information (needed for plotting)
season_info <- combined_data$season[match(rownames(community_complete),
                                              combined_data$sampling_points)]
site_info <- combined_data$site[match(rownames(community_complete),
                                      combined_data$sampling_points)]
environment_info <- combined_data$environment[match(rownames(community_complete),
                                                    combined_data$sampling_points)]

# Add labels and season information
site_scores$site_label <- rownames(site_scores)
site_scores$season <- season_info
site_scores$site <- site_info
site_scores$environment <- environment_info
species_scores$species_label <- rownames(species_scores)
env_scores$env_label <- rownames(env_scores)

# Convert environmental variable names for better visualization
env_scores$env_label <- case_when(
  env_scores$env_label == "salinity" ~ "Salinity",
  env_scores$env_label == "temperature" ~ "Temperature", 
  env_scores$env_label == "pH" ~ "pH",
  env_scores$env_label == "dissolved_o2" ~ "Dissolved O₂",
  env_scores$env_label == "sediment_particle_size" ~ "Sediment particle size",
  env_scores$env_label == "suspended_solids" ~ "Suspended solids",
  TRUE ~ env_scores$env_label
)

# Get dbRDA axis labels with variance explained
dbrda_summary <- summary(dbrda_result)
axis1_var <- round(dbrda_summary$cont$importance[2,1] * 100, 1)
axis2_var <- round(dbrda_summary$cont$importance[2,2] * 100, 1)

dbrda_plot <- ggplot() +
  geom_segment(data = env_scores, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               color = "darkgreen", arrow = arrow(length = unit(1, "cm")),
               linewidth = 1.5, alpha = 0.7) +
  geom_text(data = env_scores, aes(x = dbRDA1, y = dbRDA2, label = env_label), 
            color = "darkgreen", size = 4, fontface = "bold", 
            vjust = -1.0, hjust = 0.5, alpha = 0.7) +
  geom_point(data = site_scores, aes(x = dbRDA1, y = dbRDA2, 
                                     color = environment, shape = season), 
             size = 4, alpha = 0.8, stroke = 1.5) +
  #stat_ellipse(data = site_scores, aes(x = dbRDA1, y = dbRDA2, color = environment), 
  #             level = 0.95, linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_text(data = site_scores, aes(x = dbRDA1, y = dbRDA2, label = site_label),
            size = 3, vjust = -0.8, hjust = 0.5, color = "black") +
  scale_color_manual(values = c("Serinhaem" = "#18b798", "Sorojo" = "#1296ff", "Marau" = "#da4ef3"),
                     name = "River Basin") +
  scale_shape_manual(values = c("Wet" = 16, "Dry" = 17),
                     name = "Season") +
  labs(title = "dbRDA Biplot - Sites and Environmental Variables",
       subtitle = paste0("Overall dbRDA significance: p = ", round(dbrda_significance$`Pr(>F)`[1], 3)),
       x = paste0("dbRDA1 (", axis1_var, "%)"),
       y = paste0("dbRDA2 (", axis2_var, "%)")) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.box = "horizontal"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "gray50") +
  guides(color = guide_legend(override.aes = list(size = 5)),
         shape = guide_legend(override.aes = list(size = 5)))

# Save the plot
#ggsave("results/dbrda_biplot_ggplot.png", dbrda_plot, width = 12, height = 8, dpi = 1200)

################################################################################
############################### SIMPER Analysis ################################
################################################################################

# SIMPER on season
simper_season <- simper(community_complete, season_info, permutations = 9999)
simper_season_results <- summary(simper_season)

# SIMPER on river basin
simper_environment <- simper(community_complete, environment_info, permutations = 9999)
simper_environment_results <- summary(simper_environment)

# Function to extract SIMPER results to data frame for CSV export
extract_simper_to_dataframe <- function(simper_results, analysis_type) {
  all_results <- data.frame()
  
  # Loop through all comparisons in the SIMPER results
  for (comparison_name in names(simper_results)) {
    comparison_data <- simper_results[[comparison_name]]
    
    # Create data frame for this comparison
    comparison_df <- data.frame(
      Analysis_Type = analysis_type,
      Comparison = comparison_name,
      Species = rownames(comparison_data),
      Average_Dissimilarity = comparison_data$average,
      Standard_Deviation = comparison_data$sd,
      Ratio = comparison_data$ratio,
      Average_Abundance_Group1 = comparison_data$ava,
      Average_Abundance_Group2 = comparison_data$avb,
      Cumulative_Contribution = comparison_data$cumsum,
      P_value = comparison_data$p,
      stringsAsFactors = FALSE
    )
    
    # Add to overall results
    all_results <- rbind(all_results, comparison_df)
  }
  
  return(all_results)
}

# Extract SIMPER results to properly formatted data frames
simper_season_df <- extract_simper_to_dataframe(simper_season_results, "Season")
simper_environment_df <- extract_simper_to_dataframe(simper_environment_results, "River_Basin")

# Combine all SIMPER results into one comprehensive data frame
simper_all_results <- rbind(simper_season_df, simper_environment_df)

# Save SIMPER results to CSV
#write.csv(simper_season_df, "results/simper_season_results.csv", row.names = FALSE)
#write.csv(simper_environment_df, "results/simper_environment_results.csv", row.names = FALSE)
#write.csv(simper_all_results, "results/simper_all_results.csv", row.names = FALSE)

# Function to create lollipop plot from SIMPER results
create_simper_lollipop <- function(simper_results, title_text, top_n = 10) {
  # Extract the first comparison (typically the main one)
  comparison_data <- simper_results[[1]]
  
  # Get top n species
  top_species <- head(comparison_data, top_n)
  
  # Prepare data for plotting
  plot_data <- data.frame(
    species = rownames(top_species),
    contribution = top_species$average * 100,  # Convert to percentage
    p_value = top_species$p,
    stringsAsFactors = FALSE
  )
  
  # Reorder species by contribution (descending)
  plot_data$species <- factor(plot_data$species, levels = plot_data$species[order(plot_data$contribution)])
  
  # Create lollipop plot
  plot <- ggplot(plot_data, aes(x = species, y = contribution)) +
    geom_segment(aes(x = species, xend = species, y = 0, yend = contribution), 
                 color = "gray40", linewidth = 1.5) +
    geom_point(aes(fill = p_value < 0.05), size = 4, color = "black", shape = 21, stroke = 1) +
    scale_fill_manual(values = c("TRUE" = "white", "FALSE" = "black"),
                      name = "p < 0.05", labels = c("No", "Yes")) +
    geom_text(aes(label = ifelse(p_value < 0.001, "   p < 0.001", paste("   p =", round(p_value, 3)))),
              hjust = -0.1, vjust = 0.5, size = 3, fontface = "bold") +
    coord_flip() +
    labs(title = title_text,
         subtitle = "SIMPER Analysis Results",
         x = "Species", 
         y = "Average Contribution (%)") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray30"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 14, face = "italic"),
      axis.text.x = element_text(size = 14, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
    ) +
    expand_limits(y = max(plot_data$contribution) * 1.15)
  
  return(plot)
}

# Create lollipop plots for both season and river basin
simper_season_plot <- create_simper_lollipop(simper_season_results, 
                                           "Top 10 Species - Seasonal Differences", 
                                           top_n = 10)

simper_environment_plot <- create_simper_lollipop(simper_environment_results, 
                                                 "Top 10 Species - River Basin Differences", 
                                                 top_n = 10)

# Combine plots side by side
simper_combined_plot <- plot_grid(
  simper_season_plot + theme(legend.position = "none"),
  simper_environment_plot + theme(legend.position = "none"),
  ncol = 2,
  align = "hv",
  labels = c("A", "B"),
  label_size = 16,
  label_fontface = "bold"
)

# Extract legend from one of the plots
simper_legend <- get_legend(
  simper_season_plot + 
    theme(legend.position = "bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10))
)

# Combine plots with legend at bottom
simper_final_panel <- plot_grid(
  simper_combined_plot,
  simper_legend,
  nrow = 2,
  rel_heights = c(1, 0.1)
)

# Save the individual plots and combined panel
#ggsave("results/simper_season_lollipop.png", simper_season_plot, width = 10, height = 8, dpi = 1200)
#ggsave("results/simper_environment_lollipop.png", simper_environment_plot, width = 10, height = 8, dpi = 1200)

#ggsave("results/simper_combined_panel.png", simper_final_panel, width = 16, height = 10, dpi = 1200)
