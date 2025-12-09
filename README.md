# Microgastropod assemblages within the little-studied Camamu Bay, Brazil

This repository contains the data analysis and code supporting the paper:  
**“Investigating spatial and seasonal differences in Microgastropod assemblages within the little-studied Camamu Bay, Brazil: a potential bioindicator for remote tropical areas?”**

This study investigates the spatial and seasonal variation in the assemblages of benthic **microgastropods (snails <5–10 mm in length)** within a poorly-studied tropical bay. As these organisms are abundant, diverse, and relatively easy to sample, with well-established taxonomy, they may prove a highly suitable group for ecological studies in such areas.

Sediment samples were collected from three river basins in Camamu Bay, Brazil, during wet and dry seasons. A total of 132 microgastropods species were recorded, demonstrating high diversity.

Statistical analysis revealed significant seasonal and spatial differences in diversity and assemblage composition, which correlated with environmental gradients. The results suggest that microgastropods are a suitable component of the biota for ecological and applied studies in marine sediments, particularly in remote tropical locations where full macrofaunal analysis may be challenging. Whilst further testing across impact gradients is needed, this approach offers a practical solution for ecological investigations in regions with limited taxonomic expertise and sampling capabilities, and highlights microgastropods as a useful indicator taxon for impact biomonitoring in tropical marine sediments.

It explores how environmental variables influence assemblage structure, diversity, and composition in benthic ecosystems of remote tropical estuaries.

---

## Repository Structure

```text
├── data/
│   ├── table_S3_environmental_variables.pdf    # Environmental variables (pH, salinity, DO, etc.)
│   ├── table_S2_microgastropods_community_abundance_data.pdf     # Species abundance matrix
│
├── results/
│   ├── statistical_analysis                    # Statistical summaries, ANOVA, PERMANOVA, VIF, GAM outputs
│   ├── plots                                   # Visualization panels (diversity, environmental, NMDS, GAM)
│
├── microgastropodes_env_analyses_eduhcgalvao.R # Main R analysis script
│
└── README.md                                   # Repository documentation
```

---

## Workflow Overview

The main script **`microgastropodes_env_analyses_eduhcgalvao.R`** performs the following analytical workflow:

1. **Data Import and Cleaning**  
   Loads and merges environmental and community datasets for 15 sampling sites.

2. **Environmental Analyses**  
   Tests spatial and seasonal differences in environmental variables such as:
   - pH, temperature, salinity, dissolved oxygen
   - Suspended solids and sediment particle size

3. **Diversity Analyses**  
   Calculates and compares:
   - Shannon–Wiener diversity
   - Number of species
   - Margalef richness
   - Pielou’s evenness  
   across seasons and river basins.

4. **Community Composition**  
   - Tests assemblage differences using **PERMANOVA** (Bray-Curtis dissimilarity).  
   - Visualizes results via **NMDS** ordination with 95% confidence ellipses.

5. **Environmental–Diversity Relationships**  
   - Fits **Generalized Additive Models (GAMs)** to explore how environmental gradients affect diversity and richness.  
   - Extracts deviance explained and significance values for each variable.

6. **Assemblage–Environment Relationships**  
   - Applies **distance-based redundancy analysis (dbRDA)** to assess environmental influence on community composition after multicollinearity filtering (VIF < 7).

7. **Species Contribution (SIMPER)**  
   - Identifies top 10 taxa contributing most to spatial and seasonal dissimilarities in community composition.

---

## R Packages Used

| Package | Functionality |
|----------|----------------|
| **dplyr**, **tidyr**, **reshape2** | Data cleaning and transformation |
| **ggplot2**, **cowplot** | Visualization and figure arrangement |
| **car** | Variance Inflation Factor (VIF) analysis and statistical tests |
| **vegan** | Ecological and multivariate analyses (PERMANOVA, NMDS, dbRDA, SIMPER, diversity) |
| **mgcv** | Generalized Additive Models (GAMs) |

---

## Reproducibility

All analyses were conducted in **R version 4.3.3**.

To reproduce the workflow:

```r
# Set working directory to repository root
setwd("path/to/Microgastropodes_as_environmental_indicators")

# Run analysis
source("microgastropodes_env_analyses_eduhcgalvao.R")
