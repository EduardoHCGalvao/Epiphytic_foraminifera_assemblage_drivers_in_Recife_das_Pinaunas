# R Script: epiphytic_foram_env_analyses.R
# Author: Eduardo HC Galvao - eduardohcgalvao@gmail.com
# Version: 1.0
# Date: 08/12/2025 (dd/mm/yyyy)
# Description: Analyses to investigate interannual and seasonal patterns in
# epiphytic foraminifera assemblages.

# The workflow includes:
# 1) Data import and cleaning
# 2) Environmental data analyses 
# 3) Assemblage analyses (composition, diversity indices and structure)
# 4) Functional groups analyses
# 5) Models to investigate relationships between environmental variables with
# diversity indices and functional groups.
# 6) SIMPER: Investigate which species contribute the most for the interannual
# and seasonal differences in community composition.

# Input data:
# 1) Assemblage data: "assemblage_data_foraminifera.csv"
# 2) Environmental variables: "environmental_variables_foraminifera.csv"

# If one the following packages are not installed run:
# install.packages("package") #nolint
library(ggplot2)     # For plotting
library(cowplot)     # For combining plots
library(ggsci)       # For color schemes
library(dplyr)       # For data manipulation
library(tidyr)       # For data reshaping
library(car)         # For additional statistical tests
library(FSA)         # For additional statistical tests
library(multcomp)    # For post-hoc tests
library(vegan)       # For diversity and statistical analysis
library(mgcv)        # For GAM models

set.seed(123) # For reproducibility

###############################################################################
############################ Load and prepare data ############################
###############################################################################

# Set working directory
setwd("C:/Users/HOUSE/Desktop/lab_geemco_eduardo/projects/PHYTAL")

# Load data
assemblage_data <- read.csv("data/assemblage_data_foraminifera.csv")
env_data <- read.csv("data/environmental_variables_foraminifera.csv")

# Convert year and functional_group to factors
assemblage_data$year <- as.factor(assemblage_data$year)
assemblage_data$functional_group <- as.factor(assemblage_data$functional_group)

# Convert to long format for monthly analysis
assemblage_data_long <- assemblage_data %>%
  pivot_longer(cols = jan:dec, names_to = "month", values_to = "abundance")

# Create monthly community matrix
community_matrix <- assemblage_data %>%
  pivot_longer(cols = jan:dec, names_to = "month", values_to = "abundance") %>%
  mutate(sample_id = paste(year, month, sep = "_")) %>%
  group_by(sample_id, year, month, taxon) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = taxon, values_from = abundance, values_fill = 0) %>%
  dplyr::select(sample_id, year, month, dplyr::everything())

# Add seasonal categories
# Wet: March–August (mar, apr, may, jun, jul, aug)
# Dry: September–February (sep, oct, nov, dec, jan, feb)
community_matrix <- community_matrix %>%
  mutate(season = case_when(
    month %in% c("mar", "apr", "may", "jun", "jul", "aug") ~ "Wet",
    month %in% c("sep", "oct", "nov", "dec", "jan", "feb") ~ "Dry",
    TRUE ~ NA_character_
  )) %>%
  mutate(season = as.factor(season))

# Transform environmental data to long format
env_long <- env_data %>%
  gather(key = "year", value = "value", -season, -abiotic_variable) %>%
  mutate(year = as.numeric(gsub("X", "", year))) %>%
  mutate(abiotic_variable = case_when(
    abiotic_variable == "Ammonia" ~ "NH4+",
    abiotic_variable == "Nitrate" ~ "NO3-",
    abiotic_variable == "Phosphate" ~ "PO4^3-",
    abiotic_variable == "Dissolved Oxygen" ~ "Dissolved O2",
    abiotic_variable == "Suspended Solids" ~ "S. Solids",
    abiotic_variable == "Temperature" ~ "SST (°C)",
    TRUE ~ abiotic_variable
  ))

# Create factor for season with proper order
env_long$season <- factor(env_long$season, levels = c("Dry", "Wet"))

# Create color palette for seasons matching FI index plot
season_colors <- c("Dry" = "#35e6c0ff", "Wet" = "#5348f1ff")

###############################################################################

# Help functions
format_p_value <- function(p) {
  if(is.na(p)) return(NA_character_)
  if(p < 0.001) return("p < 0.001")
  if(p < 0.01) return(paste0("p = ", format(round(p, 3), nsmall = 3)))
  return(paste0("p = ", format(round(p, 3), nsmall = 3)))
}

safe_div <- function(x, y) {
  if(y == 0 || is.na(y)) return(NA_real_)
  x / y
}

###############################################################################
######################### Environmental data analyses #########################
###############################################################################

# Summary statistics by season (Dry vs Wet)
seasonal_summary <- env_long %>%
  filter(!is.na(season)) %>%
  group_by(abiotic_variable, season) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

# Print seasonal summary
print(seasonal_summary)

# Prepare data for plotting (with year and season)
env_summary <- env_long %>%
  filter(!is.na(season)) %>%
  group_by(abiotic_variable, year, season) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

# Combined seasonal-interannual plot
plot_combined_abiotic <- ggplot(env_summary,
                                 aes(x = factor(year),
                                     y = mean_value,
                                     group = abiotic_variable)) +
  geom_line(aes(group = interaction(abiotic_variable, season), color = season),
            linewidth = 1, alpha = 0.4,) +
  geom_point(size = 5, fill = "#ffffff", aes(color = season), 
             shape = 16, stroke = 1) +
  facet_wrap(~ abiotic_variable, scales = "free_y", ncol = 2) +
  scale_color_manual(values = season_colors, na.translate = FALSE) +
  labs(
    x = "Year",
    y = "Value"
  ) +
  theme_minimal(base_size = 20) +
theme(
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 20, face = "bold"),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# Save plots and seasonal statistical summary
#ggsave("results/environmental_variables.png", plot = plot_combined_abiotic, width = 16, height = 14, dpi = 1200)
#write.csv(seasonal_summary, "results/seasonal_environmental_summary.csv", row.names = FALSE)

###############################################################################
############################# Assemblage analyses #############################
###############################################################################

# Function to calculate FORAM Index
calculate_fi_index <- function(data) {
  # Calculate FI index for each sample using: FI = (10 * Ps) + Po + (2 * Ph)
  fi_data <- data %>%
    group_by(year, month) %>%
    summarise(
      ph_opportunistic = sum(abundance[functional_group == "Opportunistic"], na.rm = TRUE),
      ps_symbiont = sum(abundance[functional_group == "Symbiont-bearing"], na.rm = TRUE),
      po_other = sum(abundance[functional_group == "Other small forms"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      total_forams = ph_opportunistic + ps_symbiont + po_other,
      # Calculate proportions
      prop_ph = ifelse(total_forams > 0, ph_opportunistic / total_forams, 0),
      prop_ps = ifelse(total_forams > 0, ps_symbiont / total_forams, 0),
      prop_po = ifelse(total_forams > 0, po_other / total_forams, 0),
      # If no foraminifers present, assign FI = 0 (worst environmental conditions)
      fi_index = ifelse(total_forams > 0,
                        (10 * prop_ps) + (prop_po) + (2 * prop_ph),
                        0),
      # Add seasonal categories
      season = case_when(
        month %in% c("mar", "apr", "may", "jun", "jul", "aug") ~ "Wet",
        month %in% c("sep", "oct", "nov", "dec", "jan", "feb") ~ "Dry",
        TRUE ~ NA_character_
      ),
      # Create proper month ordering
      month_ordered = factor(month, 
                            levels = c("jan", "feb", "mar", "apr", "may", "jun", 
                                      "jul", "aug", "sep", "oct", "nov", "dec"),
                            labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                      "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")),
      year = factor(year, levels = c("2000", "2010", "2022")),
      season = factor(season, levels = c("Wet", "Dry"))
    )
  return(fi_data)
}

# Function to calculate diversity indices
calc_diversity <- function(community_matrix, year_col) {
  # Extract species columns (exclude metadata and taxonomy columns)
  species_cols <- community_matrix %>%
    dplyr::select(-c(year, sample_id, any_of("month"), any_of("season"), 
                     any_of("Order"), any_of("Family"), any_of("Genus")))
  # Calculate diversity indices
  diversity_df <- data.frame(
    sample_id = community_matrix$sample_id,
    year = community_matrix[[year_col]],
    month = community_matrix$month,
    season = community_matrix$season,
    num_of_species = specnumber(species_cols), 
    shannon = diversity(species_cols, index = "shannon"), 
    margalef = (specnumber(species_cols) - 1) / log(rowSums(species_cols)),
    evenness = diversity(species_cols, index = "shannon") / log(specnumber(species_cols))
  )
  return(diversity_df)
}

# Calculate diversity
diversity_monthly <- calc_diversity(community_matrix, "year")

# Calculate FI index and add it to diversity data
fi_data_for_merge <- calculate_fi_index(assemblage_data_long) %>%
  mutate(sample_id = paste(year, month, sep = "_")) %>%
  dplyr::select(sample_id, fi_index)

# Merge FI index with diversity data
diversity_monthly <- diversity_monthly %>%
  left_join(fi_data_for_merge, by = "sample_id") %>%
  dplyr::select(-num_of_species) %>%
  relocate(fi_index, .after = season) %>%
  # Ensure FI index values are valid (replace any remaining NAs with 0)
  mutate(fi_index = ifelse(is.na(fi_index), 0, fi_index))

# Define function to compute statistics for groups
compute_diversity_stats <- function(data, group_var, group_values, diversity_list) {
  result <- data.frame()
  for (val in group_values) {
    subset_data <- data[data[[group_var]] == val, ]
    for (diversity_index in diversity_list) {
      if (diversity_index %in% colnames(subset_data)) {
        diversity_values <- subset_data[[diversity_index]]
        diversity_values <- diversity_values[!is.na(diversity_values)]
        if (length(diversity_values) > 0) {
          result <- rbind(result, data.frame(
            Group = val,
            Variable = diversity_index,
            Mean = round(mean(diversity_values), 3),
            SD = round(sd(diversity_values), 3)
          ))
        }
      }
    }
  }
  return(result)
}

# List of diversity indices
diversity_indices <- c("shannon", "fi_index", "evenness", "margalef")

# Calculate diversity stats
year_values <- unique(diversity_monthly$year)
year_stats <- compute_diversity_stats(diversity_monthly, "year", year_values, diversity_indices)
month_values <- unique(diversity_monthly$month)
month_stats <- compute_diversity_stats(diversity_monthly, "month", month_values, diversity_indices) 
season_values <- unique(diversity_monthly$season)
season_stats <- compute_diversity_stats(diversity_monthly, "season", season_values, diversity_indices)

# Save diversity stats:
#write.csv(diversity_monthly, "results/diversity_indices_monthly.csv", row.names = FALSE)
#write.csv(year_stats, "results/diversity_stats_by_year.csv", row.names = FALSE)
#write.csv(month_stats, "results/diversity_stats_by_month.csv", row.names = FALSE)
#write.csv(season_stats, "results/diversity_stats_by_season.csv", row.names = FALSE)

# Prepare data for analysis
diversity_year_comp <- diversity_monthly
species_matrix_year <- community_matrix %>%
  dplyr::select(-c(year, sample_id, month, season, 
                   any_of("Order"), any_of("Family"), any_of("Genus")))

# Define diversity variables to test
diversity_vars <- c("shannon", "fi_index", "evenness", "margalef")

######################### Temporal diversity analyses #########################

# Test normality and perform statistical tests for temporal comparison
temporal_test_results <- list()

for(var in diversity_vars) {
  data_subset <- diversity_year_comp[!is.na(diversity_year_comp[[var]]) & 
                                   is.finite(diversity_year_comp[[var]]), ]
  if(nrow(data_subset) >= 6) { # Need at least 6 samples for three groups
    # Split data by years
    groups <- split(data_subset[[var]], data_subset$year)
    group_names <- names(groups)
    # Test normality for each group
    group_normality <- list()
    for(i in seq_along(groups)) {
      if(length(groups[[i]]) >= 3) {
        shapiro_result <- shapiro.test(groups[[i]])
        group_normality[[group_names[i]]] <- shapiro_result$p.value
        cat(sprintf("Normality test for %s - %s: W = %.4f, p = %.4f %s\n",
                    var, group_names[i], shapiro_result$statistic, shapiro_result$p.value,
                    ifelse(shapiro_result$p.value < 0.05, "(non-normal)", "(normal)")))
      }
    }
    # Test homogeneity of variances (Levene's test)
    if(length(groups) >= 2 && all(sapply(groups, length) >= 2)) {
      levene_test <- car::leveneTest(data_subset[[var]], data_subset$year)
      cat(sprintf("Levene's test for %s: F = %.4f, p = %.4f %s\n",
                  var, levene_test$`F value`[1], levene_test$`Pr(>F)`[1],
                  ifelse(levene_test$`Pr(>F)`[1] < 0.05, "(heterogeneous)", "(homogeneous)")))
      # Choose appropriate test
      most_normal <- mean(unlist(group_normality), na.rm = TRUE) > 0.05
      homogeneous_variances <- levene_test$`Pr(>F)`[1] > 0.05
      if(most_normal && homogeneous_variances) {
        # One-way ANOVA
        anova_model <- aov(data_subset[[var]] ~ data_subset$year)
        anova_summary <- summary(anova_model)
        f_stat <- anova_summary[[1]][["F value"]][1]
        p_value <- anova_summary[[1]][["Pr(>F)"]][1]
        cat(sprintf("One-way ANOVA for %s: F = %.4f, p = %.4f\n", var, f_stat, p_value))
        # Post-hoc test (Tukey HSD)
        if(p_value < 0.05) {
          tukey_result <- TukeyHSD(anova_model)
          posthoc_summary <- tukey_result$`data_subset$year`
          cat("Tukey HSD post-hoc test results:\n")
          print(posthoc_summary)
        } else {
          posthoc_summary <- NULL
        }
        temporal_test_results[[var]] <- list(
          test = "One-way ANOVA",
          p_value = p_value,
          statistic = f_stat,
          posthoc = posthoc_summary)
      } else {
        # Kruskal-Wallis test
        kw_test <- kruskal.test(data_subset[[var]], data_subset$year)
        cat(sprintf("Kruskal-Wallis test for %s: H = %.4f, p = %.4f\n", 
                    var, kw_test$statistic, kw_test$p.value))
        # Post-hoc test (Dunn's test) if significant
        if(kw_test$p.value < 0.05) {
          dunn_result <- FSA::dunnTest(data_subset[[var]], data_subset$year, method = "bonferroni")
          posthoc_summary <- dunn_result$res
          cat("Dunn's post-hoc test results:\n")
          print(posthoc_summary)
        } else {
          posthoc_summary <- NULL
        }
        temporal_test_results[[var]] <- list(
          test = "Kruskal-Wallis",
          p_value = kw_test$p.value,
          statistic = kw_test$statistic,
          posthoc = posthoc_summary)
      }
    }
  }
}

# Function to create post-hoc letters
create_posthoc_letters <- function(test_results, diversity_data) {
  letters_list <- list()
  for(index_name in names(test_results)) {
    if(!is.null(test_results[[index_name]]$posthoc)) {
      index_data <- diversity_data %>% 
        dplyr::select(year, all_of(index_name)) %>%
        dplyr::rename(value = all_of(index_name))
      posthoc_data <- test_results[[index_name]]$posthoc
      years <- unique(index_data$year)
      comp_matrix <- matrix(FALSE, nrow = length(years), ncol = length(years))
      rownames(comp_matrix) <- years
      colnames(comp_matrix) <- years
      if("Comparison" %in% colnames(posthoc_data) && "P.adj" %in% colnames(posthoc_data)) {
        for(i in 1:nrow(posthoc_data)) {
          comparison <- posthoc_data[i, "Comparison"]
          p_adj <- posthoc_data[i, "P.adj"]
          comp_parts <- strsplit(as.character(comparison), " - ")[[1]]
          if(length(comp_parts) == 2 && p_adj < 0.05) {
            year1 <- comp_parts[1]
            year2 <- comp_parts[2]
            if(year1 %in% years && year2 %in% years) {
              comp_matrix[year1, year2] <- TRUE
              comp_matrix[year2, year1] <- TRUE
            }
          }
        }
      } else if("p adj" %in% colnames(posthoc_data)) {
        for(i in 1:nrow(posthoc_data)) {
          comparison <- rownames(posthoc_data)[i]
          p_adj <- posthoc_data[i, "p adj"]
          comp_parts <- strsplit(as.character(comparison), "-")[[1]]
          if(length(comp_parts) == 2 && p_adj < 0.05) {
            year1 <- comp_parts[1]
            year2 <- comp_parts[2]
            if(year1 %in% years && year2 %in% years) {
              comp_matrix[year1, year2] <- TRUE
              comp_matrix[year2, year1] <- TRUE
            }
          }
        }
      } else {
        cat(sprintf("Warning: Unexpected column structure for %s: %s\n", 
                    index_name, paste(colnames(posthoc_data), collapse = ", ")))
      }
      letter_assignments <- rep("a", length(years))
      names(letter_assignments) <- years
      current_letter <- 1
      for(i in 1:length(years)) {
        for(j in 1:length(years)) {
          if(i != j && comp_matrix[years[i], years[j]]) {
            if(letter_assignments[years[i]] == letter_assignments[years[j]]) {
              letter_assignments[years[j]] <- letters[current_letter + 1]
              current_letter <- current_letter + 1
            }
          }
        }
      }
      max_values <- index_data %>%
        group_by(year) %>%
        summarise(max_val = max(value, na.rm = TRUE), .groups = 'drop')
      y_positions <- max_values$max_val * 1.04
      names(y_positions) <- max_values$year
      letters_df <- data.frame(
        index = index_name,
        year = names(letter_assignments),
        letter = as.character(letter_assignments),
        y_position = y_positions[names(letter_assignments)],
        stringsAsFactors = FALSE
      )
      letters_list[[index_name]] <- letters_df
    } else {
      years <- unique(diversity_data$year)
      index_data <- diversity_data %>% 
        dplyr::select(year, all_of(index_name)) %>%
        dplyr::rename(value = all_of(index_name))
      max_values <- index_data %>%
        group_by(year) %>%
        summarise(max_val = max(value, na.rm = TRUE), .groups = 'drop')
      letters_df <- data.frame(
        index = index_name,
        year = as.character(years),
        letter = "a",
        y_position = max_values$max_val * 1.04,
        stringsAsFactors = FALSE
      ) 
      letters_list[[index_name]] <- letters_df
    }
  }
  return(do.call(rbind, letters_list))
}

# Create post-hoc letters dataframe
posthoc_letters_df <- create_posthoc_letters(temporal_test_results, diversity_year_comp)

# Ensure proper ordering of post-hoc letters
posthoc_letters_df <- posthoc_letters_df %>%
  mutate(index = factor(index, levels = c("shannon", "fi_index", "evenness", "margalef"))) %>%
  arrange(index, year)

# Create temporal diversity boxplot
diversity_temporal_long <- diversity_year_comp %>%
  dplyr::select(sample_id, year, shannon, fi_index, evenness, margalef) %>%
  pivot_longer(cols = c(shannon, fi_index, evenness, margalef), 
               names_to = "index", values_to = "value") %>%
  mutate(index = factor(index, levels = c("shannon", "fi_index", "evenness", "margalef")))

# Add p-value annotations for temporal
temporal_p_annotations <- data.frame(
  index = c("shannon", "fi_index", "evenness", "margalef"),
  p_value_text = sapply(c("shannon", "fi_index", "evenness", "margalef"), function(var) {
    if (var %in% names(temporal_test_results)) {
      p_val <- temporal_test_results[[var]]$p_value
      test_type <- temporal_test_results[[var]]$test
      if(p_val < 0.001) return("p < 0.001")
      else return(paste("p =", round(p_val, 3)))
    } else {
      return("p = NS")
    }
  }),
  test_type = sapply(c("shannon", "fi_index", "evenness", "margalef"), function(var) {
    if (var %in% names(temporal_test_results)) {
      return(temporal_test_results[[var]]$test)
    } else {
      return("NS")
    }
  }),
  panel_labels = c("A", "B", "C", "D")
) %>%
  mutate(index = factor(index, levels = c("shannon", "fi_index", "evenness", "margalef")))

# Create separate plots for each diversity index with manual y-axis limits
p_shannon_temp <- diversity_temporal_long %>% filter(index == "shannon") %>%
  ggplot(aes(x = year, y = value, fill = year)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, linewidth = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 3) +
  geom_text(data = posthoc_letters_df %>% filter(index == "shannon"), 
            aes(x = year, y = y_position, label = letter),
            inherit.aes = FALSE, size = 8, fontface = "bold", color = "black") +
  geom_text(data = temporal_p_annotations %>% filter(index == "shannon"),
            aes(x = Inf, y = Inf, label = p_value_text), 
            hjust = 1.05, vjust = 1.55, size = 8, color = "black", fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = temporal_p_annotations %>% filter(index == "shannon"),
            aes(x = -Inf, y = Inf, label = panel_labels), 
            hjust = -0.5, vjust = 1.5, size = 8, color = "black", fontface = "bold", inherit.aes = FALSE) +
  scale_fill_npg() +
  coord_cartesian(ylim = c(1.5, 3)) +
  labs(title = "Shannon-Wiener", fill = "Year") +
  theme_classic() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 18),
        legend.key.size = unit(1.5, "cm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5))

p_fi_temp <- diversity_temporal_long %>% filter(index == "fi_index") %>%
  ggplot(aes(x = year, y = value, fill = year)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, linewidth = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 3) +
  geom_text(data = posthoc_letters_df %>% filter(index == "fi_index"), 
            aes(x = year, y = y_position, label = letter),
            inherit.aes = FALSE, size = 8, fontface = "bold", color = "black") +
  geom_text(data = temporal_p_annotations %>% filter(index == "fi_index"),
            aes(x = Inf, y = Inf, label = p_value_text), 
            hjust = 1.05, vjust = 1.55, size = 8, color = "black", fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = temporal_p_annotations %>% filter(index == "fi_index"),
            aes(x = -Inf, y = Inf, label = panel_labels), 
            hjust = -0.5, vjust = 1.5, size = 8, color = "black", fontface = "bold", inherit.aes = FALSE) +
  scale_fill_npg() +
  coord_cartesian(ylim = c(5, 9)) +
  labs(title = "FORAM Index") +
  theme_classic() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5))

p_evenness_temp <- diversity_temporal_long %>% filter(index == "evenness") %>%
  ggplot(aes(x = year, y = value, fill = year)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, linewidth = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 3) +
  geom_text(data = posthoc_letters_df %>% filter(index == "evenness"), 
            aes(x = year, y = y_position, label = letter),
            inherit.aes = FALSE, size = 8, fontface = "bold", color = "black") +
  geom_text(data = temporal_p_annotations %>% filter(index == "evenness"),
            aes(x = Inf, y = Inf, label = p_value_text), 
            hjust = 1.05, vjust = 1.55, size = 8, color = "black", fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = temporal_p_annotations %>% filter(index == "evenness"),
            aes(x = -Inf, y = Inf, label = panel_labels), 
            hjust = -0.5, vjust = 1.5, size = 8, color = "black", fontface = "bold", inherit.aes = FALSE) +
  scale_fill_npg() +
  coord_cartesian(ylim = c(0.4, 0.9)) +
  labs(title = "Pielou's Evenness") +
  theme_classic() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5))

p_margalef_temp <- diversity_temporal_long %>% filter(index == "margalef") %>%
  ggplot(aes(x = year, y = value, fill = year)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, linewidth = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 3) +
  geom_text(data = posthoc_letters_df %>% filter(index == "margalef"), 
            aes(x = year, y = y_position, label = letter),
            inherit.aes = FALSE, size = 8, fontface = "bold", color = "black") +
  geom_text(data = temporal_p_annotations %>% filter(index == "margalef"),
            aes(x = Inf, y = Inf, label = p_value_text), 
            hjust = 1.05, vjust = 1.55, size = 8, color = "black", fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = temporal_p_annotations %>% filter(index == "margalef"),
            aes(x = -Inf, y = Inf, label = panel_labels), 
            hjust = -0.5, vjust = 1.5, size = 8, color = "black", fontface = "bold", inherit.aes = FALSE) +
  scale_fill_npg() +
  coord_cartesian(ylim = c(5.5, 9)) +
  labs(title = "Margalef Richness") +
  theme_classic() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5))

# Extract legend from temporal plot
legend_temporal <- get_legend(p_shannon_temp)

# Combine temporal plots with shared legend at bottom
p_diversity_temporal <- plot_grid(
  plot_grid(p_shannon_temp + theme(legend.position = "none"), p_fi_temp, p_evenness_temp, p_margalef_temp, 
            ncol = 2, align = "hv"),
  legend_temporal,
  ncol = 1, rel_heights = c(1, 0.1)
)

######################### Seasonal diversity analyses #########################

# Test normality and perform statistical tests for seasonal comparison (Wet vs Dry)
seasonal_test_results <- list()

for(var in diversity_vars) {
    data_subset <- diversity_year_comp[!is.na(diversity_year_comp[[var]]) & 
                                   is.finite(diversity_year_comp[[var]]), ]
  if(nrow(data_subset) >= 3) {
    groups <- split(data_subset[[var]], data_subset$season)
    group_names <- names(groups)
    group_normality <- list()
    for(i in seq_along(groups)) {
      if(length(groups[[i]]) >= 3) {
        shapiro_result <- shapiro.test(groups[[i]])
        group_normality[[group_names[i]]] <- shapiro_result$p.value
        cat(sprintf("    %s (%s): W = %.4f, p = %.4f %s\n", 
                    var, group_names[i], shapiro_result$statistic, shapiro_result$p.value,
                    ifelse(shapiro_result$p.value < 0.05, "(non-normal)", "(normal)")))
      }
    }
    if(length(groups) == 2 && all(sapply(groups, length) >= 2)) {
      var_test <- var.test(groups[[1]], groups[[2]])
      cat(sprintf("  Equality of variances: F = %.4f, p = %.4f %s\n", 
                  var_test$statistic, var_test$p.value,
                  ifelse(var_test$p.value < 0.05, "(unequal)", "(equal)")))
      both_normal <- all(unlist(group_normality) > 0.05)
      equal_variances <- var_test$p.value > 0.05
      if(both_normal && equal_variances) {
        # Student's t-test
        test <- t.test(groups[[1]], groups[[2]], var.equal = TRUE)
        seasonal_test_results[[var]] <- list(
          test = "Student's t-test",
          p_value = test$p.value,
          statistic = test$statistic)
      } else if(both_normal && !equal_variances) {
        # Welch's t-test
        test <- t.test(groups[[1]], groups[[2]], var.equal = FALSE)
        seasonal_test_results[[var]] <- list(
          test = "Welch's t-test",
          p_value = test$p.value,
          statistic = test$statistic)
      } else {
        # Mann-Whitney U test
        test <- wilcox.test(groups[[1]], groups[[2]])
        n1 <- length(groups[[1]])
        n2 <- length(groups[[2]])
        U <- test$statistic
        seasonal_test_results[[var]] <- list(
          test = "Mann-Whitney U",
          p_value = test$p.value,
          statistic = test$statistic)
      }
    }
  }
}

# Create seasonal diversity boxplot
diversity_seasonal_long <- diversity_year_comp %>%
  dplyr::select(sample_id, season, shannon, fi_index, evenness, margalef) %>%
  pivot_longer(cols = c(shannon, fi_index, evenness, margalef), 
               names_to = "index", values_to = "value") %>%
  mutate(index = factor(index, levels = c("shannon", "fi_index", "evenness", "margalef")))

# Add p-value annotations for seasonal
seasonal_p_annotations <- data.frame(
  index = c("shannon", "fi_index", "evenness", "margalef"),
  p_value_text = sapply(diversity_vars, function(var) {
    if (var %in% names(seasonal_test_results)) {
      p_val <- seasonal_test_results[[var]]$p_value
      if(p_val < 0.001) return("p < 0.001")
      else return(paste("p =", round(p_val, 3)))
    } else {
      return("p = NS")
    }
  }),
  panel_labels = c("E", "F", "G", "H")
) %>%
  mutate(index = factor(index, levels = c("shannon", "fi_index", "evenness", "margalef")))

# Create separate seasonal plots for each diversity index with manual y-axis limits
p_shannon_seas <- diversity_seasonal_long %>% filter(index == "shannon") %>%
  ggplot(aes(x = season, y = value, fill = season)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, linewidth = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 3) +
  geom_text(data = seasonal_p_annotations %>% filter(index == "shannon"),
            aes(x = Inf, y = Inf, label = p_value_text), 
            hjust = 1.05, vjust = 1.55, size = 8, color = "black", fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = seasonal_p_annotations %>% filter(index == "shannon"),
            aes(x = -Inf, y = Inf, label = panel_labels), 
            hjust = -0.5, vjust = 1.5, size = 8, color = "black", fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Dry" = "#35e6c0ff", "Wet" = "#5348f1ff")) +
  coord_cartesian(ylim = c(1.5, 3)) +
  labs(title = "Shannon-Wiener") +
  theme_classic() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5))

p_fi_seas <- diversity_seasonal_long %>% filter(index == "fi_index") %>%
  ggplot(aes(x = season, y = value, fill = season)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, linewidth = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 3) +
  geom_text(data = seasonal_p_annotations %>% filter(index == "fi_index"),
            aes(x = Inf, y = Inf, label = p_value_text), 
            hjust = 1.05, vjust = 1.55, size = 8, color = "black", fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = seasonal_p_annotations %>% filter(index == "fi_index"),
            aes(x = -Inf, y = Inf, label = panel_labels), 
            hjust = -0.5, vjust = 1.5, size = 8, color = "black", fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Dry" = "#35e6c0ff", "Wet" = "#5348f1ff")) +
  coord_cartesian(ylim = c(5, 9)) +
  labs(title = "FORAM Index") +
  theme_classic() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5))

p_evenness_seas <- diversity_seasonal_long %>% filter(index == "evenness") %>%
  ggplot(aes(x = season, y = value, fill = season)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, linewidth = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 3) +
  geom_text(data = seasonal_p_annotations %>% filter(index == "evenness"),
            aes(x = Inf, y = Inf, label = p_value_text), 
            hjust = 1.05, vjust = 1.55, size = 8, color = "black", fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = seasonal_p_annotations %>% filter(index == "evenness"),
            aes(x = -Inf, y = Inf, label = panel_labels), 
            hjust = -0.5, vjust = 1.5, size = 8, color = "black", fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Dry" = "#35e6c0ff", "Wet" = "#5348f1ff")) +
  coord_cartesian(ylim = c(0.4, 0.9)) +
  labs(title = "Pielou's Evenness", fill = "Season") +
  theme_classic() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 18),
        legend.key.size = unit(1.5, "cm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5))

p_margalef_seas <- diversity_seasonal_long %>% filter(index == "margalef") %>%
  ggplot(aes(x = season, y = value, fill = season)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, linewidth = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 3) +
  geom_text(data = seasonal_p_annotations %>% filter(index == "margalef"),
            aes(x = Inf, y = Inf, label = p_value_text), 
            hjust = 1.05, vjust = 1.55, size = 8, color = "black", fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = seasonal_p_annotations %>% filter(index == "margalef"),
            aes(x = -Inf, y = Inf, label = panel_labels), 
            hjust = -0.5, vjust = 1.5, size = 8, color = "black", fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Dry" = "#35e6c0ff", "Wet" = "#5348f1ff")) +
  coord_cartesian(ylim = c(5.5, 9)) +
  labs(title = "Margalef Richness") +
  theme_classic() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5))

# Extract legend from one of the seasonal plots
legend_seasonal <- get_legend(p_evenness_seas)

# Combine seasonal plots with shared legend at bottom
p_diversity_seasonal <- plot_grid(
  plot_grid(p_shannon_seas, p_fi_seas, p_evenness_seas + theme(legend.position = "none"), p_margalef_seas, 
            ncol = 2, align = "hv"),
  legend_seasonal,
  ncol = 1, rel_heights = c(1, 0.1)
)

# Combine plots - diversity boxplots on left
diversity_combined_plot <- plot_grid(p_diversity_temporal, p_diversity_seasonal, ncol = 2)

# Save plot
#ggsave("results/diversity_indices_temporal_seasonal.png", diversity_combined_plot, width = 18, height = 14, dpi = 1200)

###############################################################################
########################## Functional group analyses ##########################
###############################################################################

# Create sample_id and prepare data
prepare_functional_group_data <- function(comm_long) {
  fg_data <- comm_long %>%
    group_by(year, month, functional_group) %>%
    summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(sample_id = paste("phytal", year, month, sep = "_")) %>%
    mutate(season = case_when(
      month %in% c("mar", "apr", "may", "jun", "jul", "aug") ~ "Wet",
      month %in% c("sep", "oct", "nov", "dec", "jan", "feb") ~ "Dry",
      TRUE ~ NA_character_
    )) %>%
    mutate(season = as.factor(season),
           year = as.factor(year),
           functional_group = as.factor(functional_group))
  return(fg_data)
}

# Summary statistics by group
summary_stats_by <- function(df, group_vars, value_var = "total_abundance") {
  value_sym <- rlang::sym(value_var)
  df %>%
    group_by(!!!rlang::syms(group_vars)) %>%
    summarise(
      mean_abundance = mean(!!value_sym, na.rm = TRUE),
      sd_abundance = sd(!!value_sym, na.rm = TRUE),
      se_abundance = sd_abundance / sqrt(n()),
      n_samples = n(),
      .groups = "drop"
    )
}

# Run temporal test for one functional group (keeps same logic as original)
run_temporal_test <- function(subdf, time_var = "year", value_var = "total_abundance", min_samples = 6) {
  res <- list()
  res$valid <- FALSE
  if(nrow(subdf) < min_samples) return(res)
  time_sym <- rlang::sym(time_var)
  val_sym <- rlang::sym(value_var)
  # ANOVA residual normality and Levene's test
  anova_model <- aov(formula = stats::as.formula(paste(value_var, "~", time_var)), data = subdf)
  sh_test <- tryCatch(shapiro.test(residuals(anova_model)), error = function(e) list(p.value = NA, statistic = NA))
  levene <- tryCatch(car::leveneTest(formula = stats::as.formula(paste(value_var, "~", time_var)), data = subdf),
                    error = function(e) NULL)
  normal_resid <- !is.na(sh_test$p.value) && sh_test$p.value > 0.05
  homogeneous_var <- !is.null(levene) && levene$`Pr(>F)`[1] > 0.05
  # Use parametric if normal residuals and homogeneous variances
  if(normal_resid && homogeneous_var) {
    anova_sum <- summary(anova_model)[[1]]
    f_stat <- anova_sum[["F value"]][1]
    p_val <- anova_sum[["Pr(>F)"]][1]
    ss_total <- sum(anova_sum[["Sum Sq"]], na.rm = TRUE)
    ss_between <- anova_sum[["Sum Sq"]][1]
    eta2 <- safe_div(ss_between, ss_total)
    posthoc <- if(!is.na(p_val) && p_val < 0.05) {
      TukeyHSD(anova_model)[[time_var]]
    } else NULL
    res <- list(
      test = "One-way ANOVA",
      statistic = f_stat,
      p_value = p_val,
      eta_squared = eta2,
      effect_size = ifelse(is.na(eta2), NA_character_,
                           ifelse(eta2 < 0.01, "small",
                                  ifelse(eta2 < 0.06, "medium", "large"))),
      posthoc = posthoc,
      shapiro = sh_test,
      levene = levene,
      valid = TRUE
    )
  } else {
    # Non-parametric Kruskal-Wallis
    kw <- kruskal.test(formula = stats::as.formula(paste(value_var, "~", time_var)), data = subdf)
    n <- length(subdf[[value_var]])
    k <- length(unique(subdf[[time_var]]))
    eps2 <- safe_div((kw$statistic - k + 1), (n - k))
    posthoc <- if(!is.na(kw$p.value) && kw$p.value < 0.05) {
      FSA::dunnTest(subdf[[value_var]], subdf[[time_var]], method = "bonferroni")$res
    } else NULL
    res <- list(
      test = "Kruskal-Wallis",
      statistic = kw$statistic,
      p_value = kw$p.value,
      epsilon_squared = eps2,
      effect_size = ifelse(is.na(eps2), NA_character_,
                           ifelse(eps2 < 0.01, "small",
                                  ifelse(eps2 < 0.06, "medium", "large"))),
      posthoc = posthoc,
      shapiro = sh_test,
      levene = levene,
      valid = TRUE
    )
  }
  return(res)
}

# Run tests across functional groups for a given analysis type (temporal vs relative vs seasonal)
run_group_tests <- function(df, grouping_var = "year", value_var = "total_abundance", min_samples = 6) {
  fgs <- unique(df$functional_group)
  tests <- list()
  for(fg in fgs) {
    subdf <- df %>% filter(functional_group == fg)
    test_res <- run_temporal_test(subdf, time_var = grouping_var, value_var = value_var, min_samples = min_samples)
    if(test_res$valid) tests[[as.character(fg)]] <- test_res
  }
  return(tests)
}

# Create post-hoc letters (compatible with Tukey and Dunn outputs as in your original code)
create_posthoc_letters <- function(test_results, data_df, value_col = "total_abundance", time_col = "year") {
  if(length(test_results) == 0) return(data.frame())
  letters_list <- list()
  for(fg in names(test_results)) {
    tinfo <- test_results[[fg]]
    # prepare the year levels and summary
    fg_data <- data_df %>% filter(functional_group == fg)
    if(nrow(fg_data) == 0) next
    years <- sort(unique(fg_data[[time_col]]))
    # default assign all 'a'
    letters_assigned <- setNames(rep("a", length(years)), years)
    posthoc <- tinfo$posthoc
    if(!is.null(posthoc)) {
      # Handle Tukey (matrix with "p adj") or Dunn (data.frame with Comparison and P.adj)
      if(is.matrix(posthoc) && "p adj" %in% colnames(posthoc)) {
        for(i in seq_len(nrow(posthoc))) {
          comp_name <- rownames(posthoc)[i]
          p_adj <- posthoc[i, "p adj"]
          parts <- strsplit(comp_name, "-")[[1]]
          if(length(parts) == 2 && p_adj < 0.05) {
            a <- parts[1]; b <- parts[2]
            if(letters_assigned[a] == letters_assigned[b]) {
              # find next available letter not used by a
              used <- unique(unlist(strsplit(letters_assigned, "")))
              next_letter <- setdiff(letters, used)[1]
              if(is.na(next_letter)) next_letter <- sample(letters, 1)
              letters_assigned[b] <- next_letter
            }
          }
        }
      } else if(is.data.frame(posthoc) && all(c("Comparison", "P.adj") %in% colnames(posthoc))) {
        for(i in seq_len(nrow(posthoc))) {
          comp <- posthoc$Comparison[i] %>% as.character()
          p_adj <- posthoc$P.adj[i]
          parts <- strsplit(comp, " - ")[[1]]
          if(length(parts) == 2 && p_adj < 0.05) {
            a <- parts[1]; b <- parts[2]
            if(letters_assigned[a] == letters_assigned[b]) {
              used <- unique(unlist(strsplit(letters_assigned, "")))
              next_letter <- setdiff(letters, used)[1]
              if(is.na(next_letter)) next_letter <- sample(letters, 1)
              letters_assigned[b] <- next_letter
            }
          }
        }
      }
    }
    # Compute plotting y positions (mean + se)*1.15
    fg_plot_summary <- fg_data %>%
      group_by_at(time_col) %>%
      summarise(mean_abundance = mean(.data[[value_col]], na.rm = TRUE),
                se_abundance = sd(.data[[value_col]], na.rm = TRUE)/sqrt(n()),
                .groups = "drop") %>%
      mutate(y_position = (mean_abundance + se_abundance) * 1.25)
    # Build letters df
    letters_df <- fg_plot_summary %>%
      mutate(functional_group = fg,
             letter = letters_assigned[as.character(.data[[time_col]])]) %>%
      dplyr::select(functional_group, all_of(time_col), letter, y_position) %>%
      rename(year = all_of(time_col))
    
    letters_list[[fg]] <- letters_df
  }
  if(length(letters_list) == 0) return(data.frame())
  out <- bind_rows(letters_list)
  out$year <- as.factor(out$year)
  return(out)
}

# Build stat annotation dataframe for plots (p-values + significance)
build_stat_annotations <- function(test_results, groups_in_plot) {
  df <- data.frame(functional_group = groups_in_plot, stringsAsFactors = FALSE)
  df$p_value <- sapply(df$functional_group, function(fg) {
    if(fg %in% names(test_results)) test_results[[fg]]$p_value else NA_real_
  })
  df$test_type <- sapply(df$functional_group, function(fg) {
    if(fg %in% names(test_results)) test_results[[fg]]$test else "No test"
  })
  df$p_text <- sapply(df$p_value, function(p) ifelse(is.na(p), "insufficient data", format_p_value(p)))
  df$significance <- sapply(df$p_value, function(p) {
    if(is.na(p)) return("--")
    if(p < 0.001) return("***")
    else if(p < 0.01) return("**")
    else if(p < 0.05) return("*")
    else return("ns")
  })
  return(df)
}

functional_group_data <- prepare_functional_group_data(assemblage_data_long)

# Summary statistics (temporal)
plot_data <- functional_group_data %>%
  group_by(functional_group, year) %>%
  summarise(
    mean_abundance = mean(total_abundance, na.rm = TRUE),
    sd_abundance = sd(total_abundance, na.rm = TRUE),
    se_abundance = sd_abundance / sqrt(n()),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(functional_group, year)

functional_test_results <- run_group_tests(functional_group_data, grouping_var = "year", value_var = "total_abundance", min_samples = 6)
functional_posthoc_letters_df <- create_posthoc_letters(functional_test_results, functional_group_data, value_col = "total_abundance", time_col = "year")
all_functional_groups <- unique(plot_data$functional_group)
stat_annotations <- build_stat_annotations(functional_test_results, all_functional_groups)

# get y positions for annotations
y_positions <- plot_data %>%
  group_by(functional_group) %>%
  summarise(max_y = max(mean_abundance + se_abundance, na.rm = TRUE) * 1.3, .groups = "drop")

stat_annotations <- stat_annotations %>% left_join(y_positions, by = "functional_group")
if(!"max_y" %in% colnames(stat_annotations)) stat_annotations$max_y <- max(plot_data$mean_abundance + plot_data$se_abundance, na.rm = TRUE) * 1.3

# Plot temporal
p_functional <- ggplot(plot_data, aes(x = functional_group, y = mean_abundance, fill = year)) +
  geom_col(position = position_dodge(width = 0.9), alpha = 0.8, color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean_abundance - se_abundance, ymax = mean_abundance + se_abundance),
                position = position_dodge(width = 0.9), width = 0.2, linewidth = 0.8) +
  geom_text(data = stat_annotations, aes(x = functional_group, y = max_y, label = p_text),
            inherit.aes = FALSE, hjust = 0.5, vjust = -0.2, size = 8, fontface = "bold") +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1.5, size = 8, fontface = "bold") +
  scale_fill_npg(name = "Year") +
  labs(y = "Mean Abundance (± SE)") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 26, face = "bold"),
    axis.text = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 24)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 5)))

# Save temporal summary and test summary
functional_summary <- functional_group_data %>%
  group_by(functional_group, year) %>%
  summarise(mean_abundance = mean(total_abundance, na.rm = TRUE),
            sd_abundance = sd(total_abundance, na.rm = TRUE),
            se_abundance = sd_abundance / sqrt(n()),
            n_samples = n(), .groups = "drop")

# Build functional_test_summary similar to original
if(length(functional_test_results) > 0) {
  functional_test_summary <- data.frame(
    Functional_Group = names(functional_test_results),
    Test_Type = sapply(functional_test_results, function(x) if(!is.null(x$test)) x$test else NA_character_),
    Statistic = sapply(functional_test_results, function(x) if(!is.null(x$statistic)) round(x$statistic, 3) else NA_real_),
    P_value = sapply(functional_test_results, function(x) if(!is.null(x$p_value)) {
      if(x$p_value < 0.001) format(x$p_value, scientific = TRUE, digits = 3) else round(x$p_value, 6)
    } else NA_character_),
    Significance = sapply(functional_test_results, function(x) {
      p <- x$p_value
      if(is.null(p) || is.na(p)) return("ns")
      if(p < 0.001) return("***")
      else if(p < 0.01) return("**")
      else if(p < 0.05) return("*")
      else return("ns")
    }),
    stringsAsFactors = FALSE
  )
} else {
  functional_test_summary <- data.frame()
}

# Relative frequencies by year and their tests
functional_relative_freq <- functional_group_data %>%
  group_by(year, month) %>%
  mutate(total_sample_abundance = sum(total_abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(relative_frequency = total_abundance / total_sample_abundance * 100) %>%
  dplyr::select(year, month, functional_group, total_abundance, total_sample_abundance, relative_frequency)

functional_relative_summary <- functional_relative_freq %>%
  group_by(year, functional_group) %>%
  summarise(
    mean_relative_freq = round(mean(relative_frequency, na.rm = TRUE), 3),
    se_relative_freq = round(sd(relative_frequency, na.rm = TRUE) / sqrt(n()), 3),
    median_relative_freq = round(median(relative_frequency, na.rm = TRUE), 3),
    n_samples = n(),
    .groups = "drop"
  )

functional_relative_test_results <- run_group_tests(functional_relative_freq, grouping_var = "year", value_var = "relative_frequency", min_samples = 6)
functional_relative_posthoc_letters_df <- create_posthoc_letters(functional_relative_test_results, functional_relative_freq, value_col = "relative_frequency", time_col = "year")

# Annotations for relative plot
all_functional_groups_rel <- unique(functional_relative_summary$functional_group)
p_value_annotations <- data.frame(
  functional_group = all_functional_groups_rel,
  p_value_text = sapply(all_functional_groups_rel, function(fg) {
    if(fg %in% names(functional_relative_test_results)) {
      p_val <- functional_relative_test_results[[fg]]$p_value
      format_p_value(p_val)
    } else {
      "insufficient data"
    }
  }),
  stringsAsFactors = FALSE
)

# compute y positions
p_value_annotations$y_pos <- sapply(p_value_annotations$functional_group, function(fg) {
  rows <- functional_relative_summary[functional_relative_summary$functional_group == fg, ]
  if(nrow(rows) == 0) return(50)
  max(rows$mean_relative_freq + rows$se_relative_freq, na.rm = TRUE) * 1.4
})

# Plot relative
p_functional_relative <- ggplot(functional_relative_summary, aes(x = functional_group, y = mean_relative_freq, fill = year)) +
  geom_col(position = position_dodge(width = 0.9), alpha = 0.8, color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = pmax(0, mean_relative_freq - se_relative_freq),
                    ymax = mean_relative_freq + se_relative_freq),
                position = position_dodge(width = 0.9), width = 0.2, linewidth = 0.8) +
  geom_text(data = p_value_annotations, aes(x = functional_group, y = y_pos, label = p_value_text),
            inherit.aes = FALSE, hjust = 0.5, vjust = -0.2, size = 8, fontface = "bold") +
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1.5, size = 8, fontface = "bold") +
  scale_fill_npg(name = "Year") +
  labs(y = "Relative Frequency (%) (± SE)") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 26, face = "bold"),
    axis.text = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 24)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Temporal seasonal tests (two-level factor comparisons)
functional_seasonal_test_results <- list()
for(fg in unique(functional_group_data$functional_group)) {
  sub <- functional_group_data %>% filter(functional_group == fg)
  if(nrow(sub) >= 3 && length(unique(na.omit(sub$season))) >= 2) {
    # reuse run_temporal_test but with grouping_var = "season"
    res <- run_temporal_test(sub, time_var = "season", value_var = "total_abundance", min_samples = 3)
    if(res$valid) functional_seasonal_test_results[[fg]] <- res
  }
}

# Plot data seasonal
plot_data_seasonal <- functional_group_data %>%
  group_by(functional_group, season) %>%
  summarise(mean_abundance = mean(total_abundance, na.rm = TRUE),
            sd_abundance = sd(total_abundance, na.rm = TRUE),
            se_abundance = sd_abundance / sqrt(n()),
            n_samples = n(), .groups = "drop") %>%
  arrange(functional_group, season)

stat_annotations_seasonal <- data.frame(
  functional_group = names(functional_seasonal_test_results),
  p_value = sapply(functional_seasonal_test_results, function(x) x$p_value),
  test_type = sapply(functional_seasonal_test_results, function(x) x$test),
  p_text = sapply(functional_seasonal_test_results, function(x) format_p_value(x$p_value)),
  significance = sapply(functional_seasonal_test_results, function(x) {
    p <- x$p_value
    if(is.null(p) || is.na(p)) return("ns")
    if(p < 0.001) return("***") else if(p < 0.01) return("**") else if(p < 0.05) return("*") else return("ns")
  }),
  stringsAsFactors = FALSE
)

# y position
y_positions_seasonal <- plot_data_seasonal %>%
  group_by(functional_group) %>%
  summarise(max_y = max(mean_abundance + se_abundance, na.rm = TRUE) * 1.1, .groups = "drop")

stat_annotations_seasonal <- stat_annotations_seasonal %>% left_join(y_positions_seasonal, by = "functional_group")
if(!"max_y" %in% colnames(stat_annotations_seasonal)) stat_annotations_seasonal$max_y <- max(plot_data_seasonal$mean_abundance + plot_data_seasonal$se_abundance, na.rm = TRUE) * 1.1

p_functional_seasonal <- ggplot(plot_data_seasonal, aes(x = functional_group, y = mean_abundance, fill = season)) +
  geom_col(position = position_dodge(width = 0.9), alpha = 0.8, color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean_abundance - se_abundance, ymax = mean_abundance + se_abundance),
                position = position_dodge(width = 0.9), width = 0.2, linewidth = 0.8) +
  geom_text(data = stat_annotations_seasonal, aes(x = functional_group, y = max_y, label = p_text),
            inherit.aes = FALSE, hjust = 0.5, vjust = -0.2, size = 8, fontface = "bold") +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1.5, size = 8, fontface = "bold") +
  scale_fill_manual(values = c("Dry" = "#35e6c0ff", "Wet" = "#5348f1ff"), name = "Season") +
  labs(y = "Mean Abundance (± SE)") +
  theme_classic() +
  theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 26, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 24)) +
  guides(fill = guide_legend(override.aes = list(size = 5)))

# Seasonal relative frequency and tests
functional_relative_freq_seasonal <- functional_group_data %>%
  group_by(season, month) %>%
  mutate(total_sample_abundance = sum(total_abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(relative_frequency = total_abundance / total_sample_abundance * 100) %>%
  dplyr::select(season, month, functional_group, total_abundance, total_sample_abundance, relative_frequency)

functional_relative_summary_seasonal <- functional_relative_freq_seasonal %>%
  group_by(season, functional_group) %>%
  summarise(mean_relative_freq = round(mean(relative_frequency, na.rm = TRUE), 3),
            se_relative_freq = round(sd(relative_frequency, na.rm = TRUE) / sqrt(n()), 3),
            median_relative_freq = round(median(relative_frequency, na.rm = TRUE), 3),
            n_samples = n(), .groups = "drop")

# Run tests for seasonal relative frequencies (two-group comparisons)
functional_relative_seasonal_test_results <- list()
for(fg in unique(functional_relative_freq_seasonal$functional_group)) {
  fg_rel_data <- functional_relative_freq_seasonal %>% filter(functional_group == fg)
  if(nrow(fg_rel_data) >= 6 && length(unique(na.omit(fg_rel_data$season))) == 2) {
    # split into seasons
    groups <- split(fg_rel_data$relative_frequency, fg_rel_data$season)
    # determine test
    normal_tests <- sapply(groups, function(x) if(length(x) >= 3) tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_) else NA_real_)
    if(length(groups) == 2 && all(sapply(groups, length) >= 2)) {
      var_test <- tryCatch(var.test(groups[[1]], groups[[2]]), error = function(e) NULL)
      both_normal <- all(normal_tests > 0.05, na.rm = TRUE)
      equal_variances <- !is.null(var_test) && var_test$p.value > 0.05
      if(both_normal && equal_variances) {
        test_res <- t.test(groups[[1]], groups[[2]], var.equal = TRUE)
        pooled_sd <- sqrt(((length(groups[[1]]) - 1) * var(groups[[1]]) + (length(groups[[2]]) - 1) * var(groups[[2]])) /
                            (length(groups[[1]]) + length(groups[[2]]) - 2))
        cohens_d <- safe_div(mean(groups[[1]], na.rm = TRUE) - mean(groups[[2]], na.rm = TRUE), pooled_sd)
        functional_relative_seasonal_test_results[[fg]] <- list(test = "Student's t-test", p_value = test_res$p.value, statistic = test_res$statistic, cohens_d = cohens_d)
      } else if(both_normal && !equal_variances) {
        test_res <- t.test(groups[[1]], groups[[2]], var.equal = FALSE)
        pooled_sd <- sqrt((var(groups[[1]], na.rm = TRUE) + var(groups[[2]], na.rm = TRUE)) / 2)
        cohens_d <- safe_div(mean(groups[[1]], na.rm = TRUE) - mean(groups[[2]], na.rm = TRUE), pooled_sd)
        functional_relative_seasonal_test_results[[fg]] <- list(test = "Welch's t-test", p_value = test_res$p.value, statistic = test_res$statistic, cohens_d = cohens_d)
      } else {
        test_res <- wilcox.test(groups[[1]], groups[[2]])
        n1 <- length(groups[[1]]); n2 <- length(groups[[2]])
        U <- test_res$statistic
        rank_biserial_r <- 1 - (2 * U) / (n1 * n2)
        functional_relative_seasonal_test_results[[fg]] <- list(test = "Mann-Whitney U", p_value = test_res$p.value, statistic = test_res$statistic, rank_biserial_r = rank_biserial_r)
      }
    }
  }
}

# Build annotation table for seasonal relative plot
if(length(functional_relative_seasonal_test_results) > 0) {
  p_value_annotations_seasonal <- data.frame(
    functional_group = names(functional_relative_seasonal_test_results),
    p_value_text = sapply(names(functional_relative_seasonal_test_results), function(fg) {
      p_val <- functional_relative_seasonal_test_results[[fg]]$p_value
      format_p_value(p_val)
    }),
    stringsAsFactors = FALSE
  )
  p_value_annotations_seasonal$y_pos <- sapply(p_value_annotations_seasonal$functional_group, function(fg) {
    rows <- functional_relative_summary_seasonal[functional_relative_summary_seasonal$functional_group == fg, ]
    if(nrow(rows) == 0) return(50)
    max(rows$mean_relative_freq + rows$se_relative_freq, na.rm = TRUE) * 1.1
  })
} else {
  p_value_annotations_seasonal <- data.frame(functional_group = character(0), p_value_text = character(0), y_pos = numeric(0))
}

p_functional_relative_seasonal <- ggplot(functional_relative_summary_seasonal, aes(x = functional_group, y = mean_relative_freq, fill = season)) +
  geom_col(position = position_dodge(width = 0.9), alpha = 0.8, color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = pmax(0, mean_relative_freq - se_relative_freq), ymax = mean_relative_freq + se_relative_freq),
                position = position_dodge(width = 0.9), width = 0.2, linewidth = 0.8) +
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.5, vjust = 1.5, size = 8, fontface = "bold") +
  scale_fill_manual(values = c("Dry" = "#35e6c0ff", "Wet" = "#5348f1ff"), name = "Season") +
  labs(y = "Relative Frequency (%) (± SE)") +
  theme_classic() +
  theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 26, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 24)) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

if(nrow(p_value_annotations_seasonal) > 0) {
  p_functional_relative_seasonal <- p_functional_relative_seasonal +
    geom_text(data = p_value_annotations_seasonal, aes(x = functional_group, y = y_pos, label = p_value_text),
              inherit.aes = FALSE, hjust = 0.5, vjust = -0.2, size = 8, fontface = "bold")
}

# Combine plots with legends
combined_functional_plots <- plot_grid(
  p_functional + theme(legend.position = "none"),
  p_functional_seasonal + theme(legend.position = "none"),
  p_functional_relative + theme(legend.position = "none"),
  p_functional_relative_seasonal + theme(legend.position = "none"),
  ncol = 2, nrow = 2, align = "hv"
)

legend_year <- get_legend(p_functional + theme(legend.position = "bottom"))
legend_season <- get_legend(p_functional_seasonal + theme(legend.position = "bottom"))
combined_legends <- plot_grid(legend_year, legend_season, ncol = 2)

final_functional_combined_plot <- plot_grid(
  combined_functional_plots,
  combined_legends,
  nrow = 2,
  rel_heights = c(1, 0.1)
)

# Extract and save post-hoc results
extract_posthoc_results_safe <- function(results_list, type_tag = "temporal") {
  if(length(results_list) == 0) return(data.frame())
  out <- purrr::map_df(names(results_list), function(fg) {
    item <- results_list[[fg]]
    if(is.null(item$posthoc)) return(data.frame(Function = fg, Comparison = NA, P_adj = NA, stringsAsFactors = FALSE))
    if(is.matrix(item$posthoc)) {
      m <- item$posthoc
      data.frame(Function = fg, Comparison = rownames(m), P_adj = m[,"p adj"], stringsAsFactors = FALSE)
    } else if(is.data.frame(item$posthoc)) {
      df <- item$posthoc
      df$Function <- fg
      df <- df %>% dplyr::select(Function, everything())
      return(df)
    } else {
      return(data.frame(Function = fg, Comparison = NA, P_adj = NA, stringsAsFactors = FALSE))
    }
  })
  return(out)
}

functional_temporal_posthoc <- extract_posthoc_results_safe(functional_test_results, "temporal")
functional_seasonal_posthoc <- extract_posthoc_results_safe(functional_seasonal_test_results, "seasonal")
functional_relative_posthoc <- extract_posthoc_results_safe(functional_relative_test_results, "temporal_relative")
functional_relative_seasonal_posthoc <- extract_posthoc_results_safe(functional_relative_seasonal_test_results, "seasonal_relative")

all_functional_posthoc_results <- bind_rows(
  functional_temporal_posthoc,
  functional_seasonal_posthoc,
  functional_relative_posthoc,
  functional_relative_seasonal_posthoc
)
rownames(all_functional_posthoc_results) <- NULL

# Save outputs
#ggsave("results/functional_groups_comprehensive_analysis.png", plot = final_functional_combined_plot, width = 18, height = 16, dpi = 1200)
#write.csv(all_functional_posthoc_results, "results/functional_groups_posthoc_comparison_results.csv", row.names = FALSE)

############## Community Composition Analysis (NMDS + PERMANOVA) ##############

# NMDS for phytal environment
nmds_phytal <- metaMDS(species_matrix_year, distance = "bray", k = 2, trymax = 100)

# Temporal comparison
permanova_temporal <- adonis2(species_matrix_year ~ year, 
                             data = community_matrix, 
                             permutations = 999, 
                             method = "bray")

# Seasonal comparison
permanova_seasonal <- adonis2(species_matrix_year ~ season, 
                             data = community_matrix, 
                             permutations = 999, 
                             method = "bray")

# Combined (Season + Year)
permanova_combined <- adonis2(species_matrix_year ~ season + year, 
                             data = community_matrix, 
                             permutations = 999, 
                             method = "bray")

# Check for homogeneity of dispersions
disp_temporal <- betadisper(vegdist(species_matrix_year, method = "bray"), 
                           community_matrix$year)
permutest_disp_temporal <- permutest(disp_temporal, permutations = 999)

disp_seasonal <- betadisper(vegdist(species_matrix_year, method = "bray"), 
                           community_matrix$season)
permutest_disp_seasonal <- permutest(disp_seasonal, permutations = 999)

# Extract results
permanova_temporal_r2 <- round(permanova_temporal$R2[1], 3)
permanova_temporal_p <- permanova_temporal$`Pr(>F)`[1]
permanova_seasonal_r2 <- round(permanova_seasonal$R2[1], 3)
permanova_seasonal_p <- permanova_seasonal$`Pr(>F)`[1]
permanova_combined_season_r2 <- round(permanova_combined$R2[1], 3)
permanova_combined_season_p <- permanova_combined$`Pr(>F)`[1]
permanova_combined_year_r2 <- round(permanova_combined$R2[2], 3)
permanova_combined_year_p <- permanova_combined$`Pr(>F)`[2]

disp_temporal_p <- permutest_disp_temporal$tab$`Pr(>F)`[1]
disp_seasonal_p <- permutest_disp_seasonal$tab$`Pr(>F)`[1]

# Create NMDS plot with seasons colored and years as shapes
nmds_scores <- data.frame(
  NMDS1 = nmds_phytal$points[,1],
  NMDS2 = nmds_phytal$points[,2],
  season = community_matrix$season,
  year = community_matrix$year
)

# Format p-values for combined PERMANOVA plot
permanova_combined_season_p_text <- if(!is.na(permanova_combined_season_p) && permanova_combined_season_p < 0.001) "p < 0.001" else paste("p =", round(permanova_combined_season_p, 3))

# Create caption text from subtitle
caption_text <- paste0("Stress = ", round(nmds_phytal$stress, 3), 
                      " | PERMANOVA (Seasons + Years): R² = ", permanova_combined_season_r2, ", ", permanova_combined_season_p_text)

p_nmds <- ggplot() +
  geom_point(data = nmds_scores, aes(x = NMDS1, y = NMDS2, color = year, shape = season),
             size = 4, alpha = 0.8, stroke = 1.5) +
  stat_ellipse(data = nmds_scores, aes(x = NMDS1, y = NMDS2, color = year),
               level = 0.95, linetype = "dashed", linewidth = 1, alpha = 0.8) +
  scale_color_npg(name = "Year") +
  scale_shape_manual(values = c("Wet" = 16, "Dry" = 17), name = "Season") +
  # Add caption in top-right corner
  annotate("text", x = Inf, y = Inf, label = caption_text, 
           hjust = 1.05, vjust = 1.5, size = 7.5, fontface = "bold") +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 22, face = "bold"),
    legend.box = "horizontal",
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 24)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "gray50") +
  guides(color = guide_legend(override.aes = list(size = 5)),
         shape = guide_legend(override.aes = list(size = 5)))

# Save NMDS plot
#ggsave("results/nmds_assemblage_composition.png", plot = p_nmds, width = 16, height = 12, dpi = 1200)

############################# Taxonomic Comparisons ############################

# Prepare data for taxonomic composition analysis
# Calculate total abundance per taxon across all samples to identify top 15
taxon_totals <- assemblage_data_long %>%
  group_by(taxon) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_abundance))

# Get top 15 most abundant taxa
top_15_taxa <- taxon_totals$taxon[1:15]

# Prepare seasonal data with relative abundances for phytal environment only
# First, identify all valid year-season combinations that have data
valid_samples <- assemblage_data_long %>%
  # Add season classification based on month
  mutate(season = case_when(
    month %in% c("mar", "apr", "may", "jun", "jul", "aug") ~ "Wet",
    month %in% c("sep", "oct", "nov", "dec", "jan", "feb") ~ "Dry",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(season)) %>%
  group_by(year, season) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  filter(total_abundance > 0) %>%
  dplyr::select(year, season)

seasonal_composition <- assemblage_data_long %>%
  # Add season classification based on month
  mutate(season = case_when(
    month %in% c("mar", "apr", "may", "jun", "jul", "aug") ~ "Wet",
    month %in% c("sep", "oct", "nov", "dec", "jan", "feb") ~ "Dry",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(season)) %>%
  # Add "Others" category for taxa not in top 15
  mutate(taxon_group = ifelse(taxon %in% top_15_taxa, taxon, "Others")) %>%
  # Calculate total abundance per sample (year + season)
  group_by(year, season, taxon_group) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  # Only complete for valid samples (those that have data)
  right_join(valid_samples, by = c("year", "season")) %>%
  complete(nesting(year, season), taxon_group, fill = list(abundance = 0)) %>%
  # Calculate total abundance per sample for relative abundance calculation
  group_by(year, season) %>%
  mutate(total_sample_abundance = sum(abundance, na.rm = TRUE),
         relative_abundance = (abundance / total_sample_abundance) * 100) %>%
  ungroup() %>%
  # Create proper season and year ordering
  mutate(season = factor(season, levels = c("Dry", "Wet")),
         year = factor(year, levels = c("2000", "2010", "2022")),
         year_season = interaction(year, season, sep = " - ")) %>%
  # Ensure Others comes last in the legend
  mutate(taxon_group = factor(taxon_group, 
                             levels = c(top_15_taxa, "Others")))

# Define custom colors for each taxon (you can personalize these colors)
# Create a named vector with colors for each of the top 15 taxa plus "Others"
taxon_colors <- c(
  # You can customize these colors for each species
  "#3d49f1ff",  # Color for taxon 1
  "#4DBBD5FF",  # Color for taxon 2
  "#00A087FF",  # Color for taxon 3
  "#5e87e0ff",  # Color for taxon 4
  "#F39B7FFF",  # Color for taxon 5
  "#8491B4FF",  # Color for taxon 6
  "#386157ff",  # Color for taxon 7
  "#f74444ff",  # Color for taxon 8
  "#693202ff",  # Color for taxon 9
  "#8d7355ff",  # Color for taxon 10
  "#73f562ff",    # Color for taxon 11
  "#e7bb41ff",    # Color for taxon 12
  "#3427e4ff",    # Color for taxon 13
  "#ac4141ff",    # Color for taxon 14
  "#98D8C8",    # Color for taxon 15
  "gray60"      # Color for "Others"
)

# Assign names to colors based on the taxon_group factor levels
names(taxon_colors) <- levels(seasonal_composition$taxon_group)

# Create the vertical stacked bar plot with faceting by year
plot_taxonomic_phytal <- ggplot(seasonal_composition, aes(x = season, y = relative_abundance, fill = taxon_group)) +
  geom_col(position = "stack", alpha = 0.9, color = "white", linewidth = 0.1) +
  facet_wrap(~ year, nrow = 1, labeller = labeller(year = c("2000" = "Year 2000", "2010" = "Year 2010", "2022" = "Year 2022"))) +
  scale_fill_manual(values = taxon_colors, name = "Species") +
  labs(x = "Season", y = "Relative Abundance (%)") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 22, face = "bold"),
    strip.text = element_text(size = 22, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 17, face = "bold"),
    legend.text = element_text(size = 17, face = "italic"),
    legend.key.size = unit(0.5, "cm"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.spacing.x = unit(0.25, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  guides(fill = guide_legend(ncol = 5, byrow = TRUE, override.aes = list(color = "black", linewidth = 0.3)))

# Calculate summary statistics for top taxa in phytal environment
taxonomic_summary <- seasonal_composition %>%
  group_by(taxon_group) %>%
  summarise(
    mean_rel_abundance = mean(relative_abundance, na.rm = TRUE),
    sd_rel_abundance = sd(relative_abundance, na.rm = TRUE),
    min_rel_abundance = min(relative_abundance, na.rm = TRUE),
    max_rel_abundance = max(relative_abundance, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_rel_abundance))

# Display the top 15 taxa information
for(i in 1:15) {
  total_abund <- taxon_totals$total_abundance[i]
  taxon_name <- taxon_totals$taxon[i]
  cat(sprintf("%2d. %-25s (Total: %6.0f)\n", i, taxon_name, total_abund))
}

# Save outputs
#ggsave("results/taxonomic_composition_phytal_environment.png", plot = plot_taxonomic_phytal, width = 16, height = 12, dpi = 1200)
#write.csv(taxonomic_summary, "results/taxonomic_composition_summary.csv", row.names = FALSE)
#write.csv(seasonal_composition, "results/seasonal_taxonomic_composition_data.csv", row.names = FALSE)

##############################################################################
################ Diversity and Environmental variables models ################
##############################################################################

# Prepare diversity data at year-season level to match environmental variables
diversity_for_correlation <- diversity_monthly %>%
  dplyr::select(year, month, season, shannon, fi_index, evenness, margalef) %>%
  # Convert to year-season summary (same as environmental variables structure)
  group_by(year, season) %>%
  summarise(
    shannon = mean(shannon, na.rm = TRUE),
    fi_index = mean(fi_index, na.rm = TRUE),
    evenness = mean(evenness, na.rm = TRUE),
    margalef = mean(margalef, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  gather(key = "diversity_index", value = "diversity_value", 
         shannon, fi_index, evenness, margalef) %>%
  mutate(
    year = as.numeric(as.character(year)),
    diversity_index = factor(diversity_index, 
                           levels = c("shannon", "fi_index", "evenness", "margalef"),
                           labels = c("Shannon-Wiener", "FORAM Index", 
                                    "Pielou's Evenness", "Margalef Richness"))
  )

# Prepare environmental variables for merging
# Use env_long instead of abiotic_summary which has the data in long format
env_for_correlation <- env_long %>%
  group_by(abiotic_variable, year, season) %>%
  summarise(abiotic_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(abiotic_value)) %>%
  # Ensure season is a factor with proper levels
  mutate(season = factor(season, levels = c("Dry", "Wet")))

# Merge environmental variables and diversity data
correlation_data <- diversity_for_correlation %>%
  left_join(env_for_correlation, by = c("year", "season")) %>%
  filter(!is.na(abiotic_value), !is.na(diversity_value))

# Calculate correlation statistics for each abiotic-diversity combination
if(nrow(correlation_data) > 0) {
  correlation_stats <- correlation_data %>%
    group_by(abiotic_variable, diversity_index) %>%
    summarise(
      n_points = n(),
      correlation = cor(abiotic_value, diversity_value, use = "complete.obs"),
      .groups = "drop"
    ) %>%
    mutate(
      correlation_strength = case_when(
        abs(correlation) >= 0.7 ~ "Strong",
        abs(correlation) >= 0.5 ~ "Moderate", 
        abs(correlation) >= 0.3 ~ "Weak",
        TRUE ~ "Very Weak"
      ),
      correlation_direction = ifelse(correlation > 0, "Positive", "Negative")
    )
}

# Perform linear regression for each combination
regression_results <- list()

env_vars <- unique(correlation_data$abiotic_variable)
diversity_indices <- unique(correlation_data$diversity_index)

for(env_var in env_vars) {
  for(diversity_idx in diversity_indices) {
    
    # Filter data for this combination
    combo_data <- correlation_data %>%
      filter(abiotic_variable == env_var, diversity_index == diversity_idx)
    
    if(nrow(combo_data) >= 3) {
      # Perform linear regression
      lm_model <- lm(diversity_value ~ abiotic_value, data = combo_data)
      lm_summary <- summary(lm_model)
      
      # Extract statistics
      r_squared <- lm_summary$r.squared
      p_value <- lm_summary$coefficients[2, 4]  # p-value for slope
      slope <- lm_summary$coefficients[2, 1]
      correlation <- cor(combo_data$abiotic_value, combo_data$diversity_value, use = "complete.obs")
      
      # Store results
      regression_results[[paste(env_var, diversity_idx, sep = "_")]] <- list(
        abiotic_variable = env_var,
        diversity_index = diversity_idx,
        n_points = nrow(combo_data),
        correlation = correlation,
        r_squared = r_squared,
        slope = slope,
        p_value = p_value,
        model = lm_model
      )
      
      # Print results
      cat(sprintf("\n%s vs %s:\n", env_var, diversity_idx))
      cat(sprintf("  Correlation: r = %.3f\n", correlation))
      cat(sprintf("  R-squared: %.3f\n", r_squared))
      cat(sprintf("  Slope: %.4f\n", slope))
      cat(sprintf("  P-value: %.4f %s\n", p_value, 
                  ifelse(p_value < 0.05, "(significant)", "(not significant)")))
    }
  }
}

# Create annotation data with R² and p-values for each facet
annotation_data <- data.frame()

for(env_var in env_vars) {
  for(diversity_idx in diversity_indices) {
    
    # Get the regression results for this combination
    combo_key <- paste(env_var, diversity_idx, sep = "_")
    
    if(combo_key %in% names(regression_results)) {
      result <- regression_results[[combo_key]]
      
      # Get the data range for positioning
      combo_data <- correlation_data %>%
        filter(abiotic_variable == env_var, diversity_index == diversity_idx)
      
      if(nrow(combo_data) >= 3) {
        # Calculate position for annotation (top-left of plot)
        x_range <- range(combo_data$abiotic_value, na.rm = TRUE)
        y_range <- range(combo_data$diversity_value, na.rm = TRUE)
        
        x_pos <- x_range[1] + 0.05 * (x_range[2] - x_range[1])
        y_pos <- y_range[2] - 0.05 * (y_range[2] - y_range[1])
        
        # Format the statistics text
        r_squared_text <- sprintf("R² = %.3f", result$r_squared)
        p_text <- if(result$p_value < 0.001) {
          "p < 0.001"
        } else {
          sprintf("p = %.3f", result$p_value)
        }
        
        # Combine R² and p-value
        stats_text <- paste0(r_squared_text, "\n", p_text)
        
        # Add to annotation data
        annotation_data <- rbind(annotation_data, data.frame(
          abiotic_variable = env_var,
          diversity_index = diversity_idx,
          x_pos = x_pos,
          y_pos = y_pos,
          stats_text = stats_text,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# Create summary table of significant correlations
significant_correlations <- correlation_stats %>%
  filter(abs(correlation) >= 0.3) %>%  # Only show moderate or stronger correlations
  arrange(desc(abs(correlation))) %>%
  mutate(
    correlation_rounded = round(correlation, 3),
    interpretation = paste0(correlation_strength, " ", correlation_direction, " correlation")
  )

# Save regression analysis results
regression_summary_table <- do.call(rbind, lapply(regression_results, function(x) {
  data.frame(
    Abiotic_Variable = x$abiotic_variable,
    Diversity_Index = x$diversity_index,
    N_Points = x$n_points,
    Correlation = round(x$correlation, 4),
    R_Squared = round(x$r_squared, 4),
    Slope = round(x$slope, 4),
    P_Value = round(x$p_value, 4),
    Significance = ifelse(x$p_value < 0.001, "***",
                         ifelse(x$p_value < 0.01, "**",
                               ifelse(x$p_value < 0.05, "*", "ns"))),
    stringsAsFactors = FALSE
  )
}))

# Create comprehensive regression plot with annotations
plot_env_diversity_regression <- ggplot(correlation_data, aes(x = abiotic_value, y = diversity_value)) +
  geom_point(aes(color = season), size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.75, alpha = 0.3) +
  geom_text(data = annotation_data, 
            aes(x = x_pos, y = y_pos, label = stats_text),
            hjust = -0.25, vjust = -0.5, size = 6, fontface = "bold",
            color = "black", inherit.aes = FALSE,
            bbox = list(boxstyle = "round,pad=0.3", facecolor = "white", alpha = 0.8)) +
  facet_grid(diversity_index ~ abiotic_variable, scales = "free") +
  scale_color_manual(values = c("Dry" = "#35e6c0ff", "Wet" = "#5348f1ff"), name = "Season") +
  labs(x = "Environmental Variable Value", y = "Diversity Index Value") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 22, face = "bold"),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(0.6, "lines"),
    legend.position = "bottom"
  )

# Save outputs
#ggsave("results/abiotic_diversity_regression_analysis_with_FORAM_index.png", plot = plot_env_diversity_regression, width = 18, height = 14, dpi = 1200)
#write.csv(regression_summary_table, "results/abiotic_diversity_regression_results.csv", row.names = FALSE)
#write.csv(correlation_stats, "results/abiotic_diversity_correlation_summary.csv", row.names = FALSE)

################################################################################
################################ SIMPER Analysis ###############################
################################################################################

# Function to perform SIMPER analysis and create plots
perform_simper_analysis <- function(comm_data, species_matrix, analysis_type = "temporal") {
  # Set grouping variable based on analysis type
  if (analysis_type == "temporal") {
    group_var <- "year"
    analysis_title <- "Interannual"
  } else if (analysis_type == "seasonal") {
    group_var <- "season"
    analysis_title <- "Seasonal"
  }
  # Remove species columns with zero variance (all zeros)
  species_variance <- apply(species_matrix, 2, function(x) var(x, na.rm = TRUE))
  non_zero_species <- species_matrix[, species_variance > 0, drop = FALSE]
  # Check grouping variable exists
  if(!group_var %in% colnames(comm_data)) {
    stop(paste("Group variable", group_var, "not found in community data"))
  }
  # Get grouping factor
  group_factor <- comm_data[[group_var]]
  # Call simper with explicit namespace
  simper_result <- vegan::simper(non_zero_species, group_factor, permutations = 999)
  # Extract results - SIMPER with 3+ groups creates multiple pairwise comparisons
  # We'll use the overall average across all comparisons for the plot
  simper_summary <- summary(simper_result)
  # For the plot, we'll combine results across all comparisons
  # This gives a better overview than just showing one pairwise comparison
  all_species <- c()
  all_contributions <- c()
  all_pvalues <- c()
  for(comp_name in names(simper_summary)) {
    comp_data <- simper_summary[[comp_name]]
    all_species <- c(all_species, rownames(comp_data))
    all_contributions <- c(all_contributions, comp_data$average)
    all_pvalues <- c(all_pvalues, comp_data$p)
  }
  # Create combined dataset and aggregate by species (average across comparisons)
  combined_data <- data.frame(
    species = all_species,
    contribution = all_contributions,
    p_value = all_pvalues,
    stringsAsFactors = FALSE
  )
  # Aggregate by species (taking mean contribution and minimum p-value)
  aggregated_data <- combined_data %>%
    group_by(species) %>%
    summarise(
      contribution = mean(contribution, na.rm = TRUE),
      p_value = min(p_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(contribution))
  # Take top 15 species for clarity
  n_species <- min(15, nrow(aggregated_data))
  plot_data <- aggregated_data[1:n_species, ]
  # Use base R operations to avoid potential conflicts
  plot_data <- plot_data[order(plot_data$contribution, decreasing = TRUE), ]
  plot_data$species <- factor(plot_data$species, levels = rev(plot_data$species))
  # Create title
  if (analysis_type == "seasonal") {
    title_text <- analysis_title
  } else {
    title_text <- paste("SIMPER - ", analysis_title)
  }
  # Create the plot
  plot <- ggplot(plot_data, aes(x = species, y = contribution)) +
    geom_segment(aes(x = species, xend = species, y = 0, yend = contribution), 
                 color = "gray40", linewidth = 1.5) +
    geom_point(size = 4, color = "black", shape = 21, fill = "black", stroke = 1) +
    coord_flip() +
    labs(title = title_text,
         x = "Species", 
         y = "Average Contribution (%)") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 24, hjust = 0.5, color = "gray30"),
      axis.title = element_text(size = 22, face = "bold"),
      axis.text.y = element_text(size = 20, face = "italic"),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 20, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 20)
    ) +
    expand_limits(y = max(plot_data$contribution) * 1.15)
  # Return results
  return(list(
    analysis_type = analysis_type,
    simper_result = simper_result,
    plot_data = plot_data,
    plot = plot,
    comparison = paste("All pairwise comparisons:", paste(names(simper_summary), collapse = ", ")),
    overall_dissimilarity = attr(simper_result, "average.dissimilarity")
  ))
}

# Prepare species matrix for SIMPER (same as used for NMDS)
species_matrix_simper <- community_matrix %>%
  dplyr::select(-c(year, sample_id, month, season, 
                   any_of("Order"), any_of("Family"), any_of("Genus")))

# Run SIMPER analysis for temporal comparisons
phytal_simper_temporal <- perform_simper_analysis(community_matrix, species_matrix_simper, "temporal")
# Phytal environment - seasonal
phytal_simper_seasonal <- perform_simper_analysis(community_matrix, species_matrix_simper, "seasonal")
# Extract plots
plot_phytal_simper_temporal <- phytal_simper_temporal$plot
plot_phytal_simper_seasonal <- phytal_simper_seasonal$plot

###################### Pairwise SIMPER Analysis Between Years ######################

# Function to perform pairwise SIMPER analysis between specific years
perform_pairwise_simper <- function(comm_data, species_matrix, year1, year2) {
  # Filter data for the two years
  pairwise_data <- comm_data %>%
    filter(year %in% c(year1, year2))
  # Filter species matrix accordingly
  pairwise_indices <- which(comm_data$year %in% c(year1, year2))
  pairwise_species <- species_matrix[pairwise_indices, ]
  # Remove species with zero abundance across both years
  species_sums <- colSums(pairwise_species)
  non_zero_species <- pairwise_species[, species_sums > 0]
  # Create grouping factor
  group_factor <- as.factor(pairwise_data$year)
  # Run SIMPER
  simper_result <- vegan::simper(non_zero_species, group_factor, permutations = 999)
  # Get summary
  simper_summary <- summary(simper_result)
  # Extract data for the comparison
  comparison_name <- names(simper_summary)[1]
  simper_data <- simper_summary[[comparison_name]]
  # Get overall dissimilarity from the simper result
  # Extract the overall dissimilarity properly from the summary
  overall_dissim <- round(sum(simper_data$average), 3)
  # Create summary table (all species)
  n_species <- nrow(simper_data)
  pairwise_table <- data.frame(
    species = rownames(simper_data)[1:n_species],
    avg_contribution = round(simper_data$average[1:n_species], 3),
    cumulative_contribution = round(simper_data$cumsum[1:n_species], 3),
    p_value = round(simper_data$p[1:n_species], 3),
    significance = ifelse(simper_data$p[1:n_species] < 0.001, "***",
                         ifelse(simper_data$p[1:n_species] < 0.01, "**",
                               ifelse(simper_data$p[1:n_species] < 0.05, "*", "ns")))
  )
  return(list(
    comparison = paste(year1, "vs", year2),
    overall_dissimilarity = overall_dissim,
    table = pairwise_table,
    simper_result = simper_result
  ))
}

# Perform all pairwise comparisons
pairwise_2000_2010 <- perform_pairwise_simper(community_matrix, species_matrix_simper, 2000, 2010)
pairwise_2000_2022 <- perform_pairwise_simper(community_matrix, species_matrix_simper, 2000, 2022)
pairwise_2010_2022 <- perform_pairwise_simper(community_matrix, species_matrix_simper, 2010, 2022)

#################### Create Individual Pairwise SIMPER Plots ####################

# Function to create SIMPER plot from pairwise results
create_pairwise_simper_plot <- function(pairwise_result, max_species = 10) {
  if(is.null(pairwise_result)) return(NULL)
  # Get top contributing species
  plot_data <- pairwise_result$table[1:min(max_species, nrow(pairwise_result$table)), ]
  # Prepare data for plotting
  plot_data$species <- factor(plot_data$species, levels = rev(plot_data$species))
  # Create the plot
  p <- ggplot(plot_data, aes(x = species, y = avg_contribution)) +
    geom_segment(aes(x = species, xend = species, y = 0, yend = avg_contribution), 
                 color = "gray40", linewidth = 1.5) +
    geom_point(size = 4, color = "black", fill = "black", shape = 21, stroke = 1) +
    coord_flip() +
    labs(
      title = pairwise_result$comparison,
      x = "Species",
      y = "Average Contribution (%)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 22, face = "bold"),
      axis.text.y = element_text(size = 20, face = "italic"),
      axis.text.x = element_text(size = 20, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    expand_limits(y = max(plot_data$avg_contribution) * 1.15)
  return(p)
}

# Create individual plots for each pairwise comparison
plot_2000_2010 <- create_pairwise_simper_plot(pairwise_2000_2010, max_species = 8)
plot_2000_2022 <- create_pairwise_simper_plot(pairwise_2000_2022, max_species = 8)
plot_2010_2022 <- create_pairwise_simper_plot(pairwise_2010_2022, max_species = 8)

# Remove axis labels for specific panels in 2x2 layout
# A (top-left): Remove x-axis label
# B (top-right): Remove x-axis and y-axis labels
# C (bottom-left): Keep both labels
# D (bottom-right): Remove y-axis label
plot_2000_2010 <- plot_2000_2010 + labs(y = NULL)  # A: Remove x-axis
plot_2000_2022 <- plot_2000_2022 + labs(x = NULL, y = NULL)  # B: Remove both x and y

# Create 2x2 panel with three pairwise comparisons and seasonal analysis
combined_simper_panel <- plot_grid(
  plot_2000_2010,
  plot_2000_2022,
  plot_2010_2022,
  plot_phytal_simper_seasonal,
  ncol = 2, nrow = 2,
  align = "hv",
  labels = c("A", "B", "C", "D"),
  label_size = 24, 
  label_fontface = "bold"
)

# Combine plots for phytal environment (original 1x2 layout for backward compatibility)
combined_phytal_simper <- plot_grid(plot_phytal_simper_temporal, plot_phytal_simper_seasonal,
                                   ncol = 2, align = "h",
                                   label_size = 24, label_fontface = "bold")

# Compile SIMPER results into summary tables
compile_simper_results <- function(simper_list) {
  results_df <- data.frame()
  for(result in simper_list) {
    env_name <- result$environment
    analysis_type <- result$analysis_type
    plot_data <- result$plot_data
    # Add environment and analysis type info
    temp_df <- plot_data %>%
      mutate(Analysis_Type = analysis_type,
             Rank = 1:nrow(plot_data),
             Significance = ifelse(p_value < 0.001, "***",
                                  ifelse(p_value < 0.01, "**",
                                        ifelse(p_value < 0.05, "*", "ns")))) %>%
      dplyr::select(Analysis_Type, Rank, Species = species, 
             Contribution = contribution, P_value = p_value, Significance)
    results_df <- rbind(results_df, temp_df)
  }
  return(results_df)
}

# Compile all SIMPER results
all_simper_results <- list(phytal_simper_temporal, phytal_simper_seasonal)
simper_summary_table <- compile_simper_results(all_simper_results)

# Compile pairwise results into detailed tables
compile_pairwise_results <- function(pairwise_list) {
  all_results <- data.frame()
  for(result in pairwise_list) {
    if(!is.null(result)) {
      # Ensure overall_dissimilarity is a single value
      overall_dissim_value <- as.numeric(result$overall_dissimilarity)[1]
      temp_df <- result$table %>%
        mutate(comparison = result$comparison,
               overall_dissimilarity = overall_dissim_value) %>%
        dplyr::select(comparison, overall_dissimilarity, species, avg_contribution, 
               cumulative_contribution, p_value, significance)
      all_results <- rbind(all_results, temp_df)
    }
  }
  return(all_results)
}

# Create detailed pairwise results table
pairwise_results_list <- list(pairwise_2000_2010, pairwise_2000_2022, pairwise_2010_2022)
pairwise_simper_table <- compile_pairwise_results(pairwise_results_list)

# Save results
#ggsave("results/simper_analysis_pairwise_panel_2x2.png", plot = combined_simper_panel, width = 20, height = 16, dpi = 1200)
#write.csv(simper_summary_table, "results/simper_analysis_results.csv", row.names = FALSE)
#write.csv(pairwise_simper_table, "results/pairwise_simper_dissimilarity_results.csv", row.names = FALSE)
