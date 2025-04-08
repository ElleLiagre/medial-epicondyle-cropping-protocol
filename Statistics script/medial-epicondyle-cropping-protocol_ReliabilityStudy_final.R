# ===========================================================================
# Title: A Standardized, Three-Dimensional Cropping Protocol for Analyzing 
#        the Medial Epicondyle of the Humerus - Reliability Study
# 
# Author: Elle B. K. Liagre
# Email: elle.liagre@u-bordeaux.fr
# ORCID: https://orcid.org/0000-0002-8993-3266
# Version: 1.0
# Date: 2024-12-03
# 
# Description: This R script includes the statistics used to test the repeatability 
# and reproducibility of a semi-automated, standardized 3D cropping protocol for 
# analyzing the medial epicondyle of the human humerus. The method ensures high 
# repeatability and reproducibility, facilitating the study of entheseal changes 
# and skeletal adaptations in biological anthropology.
#
#
# License: This script is licensed under the GNU General Public License v3.0.
# See https://www.gnu.org/licenses/gpl-3.0.html for more details.
#  ===========================================================================
#
# Required Packages: 
#   - dplyr (>= 1.1.4)
#   - tidyverse (>= 2.0.0)
#   - purrr (>= 1.0.2)
#   - openxlsx (>= 4.2.5.2)
#   - GUniFrac (>= 1.8)
#   - irr (>= 0.84.1)
#   - DescTools (>= 0.99.57)
#   - ggplot2 (>= 3.4.4)
#   - tidyr (>= 1.3.1)
#   - reshape2 (>= 1.4.4)
#   - RColorBrewer (>=1.1-3)
# 
# Input Data:
#   - Two CSV files containing:
#       1. Raw 3D coordinates of landmarks
#       2. Raw 3D surface areas of the cropped models
# These files can be found as supplementary information accompanying the article. 
#     
# 
# Output Data:
#   - Excel files:
#       1. Landmark data: Euclidean distances between corresponding points
#       2. Landmark data: distance-based Intraclass Correlation Coefficient (dICC) results
#       3. 3D Surface Area data: Intra-observer results (descriptive statistics, 
#          percentage error)
#       4. 3D Surface Area data: Inter-observer results (descriptive statistics, 
#          percentage error)
#       5. 3D Surface Area data: Correlation coefficients (normality test, 
#          Intraclass Correlation Coefficient (ICC) results, Lin's Concordance
#          Correlation Coefficient (CCC) results and correlation matrix)
#   - Figures:
#       1. Boxplot of Euclidean distances between corresponding points
#       2. Lin's CCC pairwise comparison matrix (heatmap)
#       3. Lin's concordance correlation plots
# 
# R Version: v4.3.2
# ===========================================================================

## Please clear environment `rm(list = ls())` and run the whole script ##


## Start of Analysis ##


# Dependencies ------------------------------------------------------------

library(dplyr)
library(tidyverse)
library(purrr)
library(openxlsx)
library(GUniFrac)
library(irr)
library(DescTools)
library(ggplot2)
library(tidyr)
library(reshape2)
library(RColorBrewer)


# Functions ---------------------------------------------------------------

# Function to calculate Euclidean distances
calculate_euclidean_distance <- function(coord1, coord2) {
  sqrt(sum((coord1 - coord2)^2))
}


# Function to calculate Mean and Standard Deviation for each error type
calculate_mean_sd <- function(error_data) {
  error_data %>%
    mutate(across(starts_with("Set"), as.numeric)) %>%
    rowwise() %>%
    mutate(
      Mean_Error = round(mean(c_across(starts_with("Set")), na.rm = TRUE), 2),
      SD_Error = round(sd(c_across(starts_with("Set")), na.rm = TRUE), 2)
    ) %>%
    ungroup()
}

# Function to append overall results of descriptive statistics
append_overall_row <- function(error_data) {
  overall_mean_error <- round(mean(error_data$Mean_Error, na.rm = TRUE), 2)
  overall_sd_error <- round(sd(error_data$Mean_Error, na.rm = TRUE), 2)
  
  overall_row <- tibble(
    Model = "Overall",
    Mean_Error = overall_mean_error,
    SD_Error = overall_sd_error
  )
  
  return(bind_rows(error_data, overall_row))
}

# Function to calculate ICC
calculate_icc <- function(data, lm_name, x_col, y_col, z_col, strata, output_df, B = 100000) {
  # Extract coordinates
  coords <- data[, c(x_col, y_col, z_col)]
  dist.mat <- as.matrix(dist(coords))
  
  # Calculate the distance-based ICC with standard error using bootstrap
  dICC_result <- dICC(dist.mat, strata)
  dICC_se_bt <- dICC.SE.bt(dist.mat, strata, B = B) 
  
  # Calculate dICC value, SE, and confidence intervals
  dICC_value <- dICC_result$ICC
  SE <- dICC_se_bt$SE
  Z <- 1.96
  
  lower_bound <- dICC_value - Z * SE
  upper_bound <- min(dICC_value + Z * SE, 1)
  
  # Round the values to 3 decimal places
  dICC_value <- round(dICC_value, 3)
  SE <- round(SE, 3)
  lower_bound <- round(lower_bound, 3)
  upper_bound <- round(upper_bound, 3)
  
  # Create a new row for the output data frame
  new_row <- data.frame(
    Inter_obs_LMx = lm_name,
    ICC = dICC_value,
    SE = SE,
    Lower_CI = lower_bound,
    Upper_CI = upper_bound,
    Iterations = B, 
    stringsAsFactors = FALSE
  )
  
  # Append the new row to the output data frame
  output_df <- rbind(output_df, new_row)
  
  return(output_df)  
}


# Function to calculate percentage error
percent_error <- function(obs1, obs2) {
  return(abs(obs1 - obs2) / ((obs1 + obs2)/2) * 100)
}

# Import data -------------------------------------------------------------

## Set the working directory to where your files are stored.
# This can be done programmatically: 
#   Right-click the folder or file, choose "Copy as path", 
#   and replace backslashes (\) with forward slashes (/).
#   Replace this path with your actual folder path:
setwd("C:/Users/Username/path/to/your/folder")

# Or, if working in RStudio, this can be done using the GUI:
#   1. Go to the "Session" menu in the top bar.
#   2. Select "Set Working Directory" and then "Choose Directory..."
#   3. Navigate to the folder you want and click "Open".
# The working directory will be set to the selected folder.

# Verify the working directory is set correctly:
getwd()


### Landmark data ###

# Import file with all landmark coordinates
coord_full <- read.csv(file.path(getwd(), "SuppInfo_4.csv"), 
                       header=TRUE, 
                       stringsAsFactors=FALSE, 
                       sep = ',')
coord_full <- coord_full %>%
  mutate(across(starts_with("LM"), ~as.numeric(gsub(",", ".", .))))

# Create intra-observer dataframe
coord_intra <- coord_full[coord_full$Observer == 1 & coord_full$Set != "Set average",]
coord_intra$Observer <- NULL

# Create inter-observer dataframe
coord_inter <- coord_full[coord_full$Observer == 2 | coord_full$Set == "Set average",]
coord_inter$Set <- NULL

# Clean up full coordinate file (remove average set)
coord_full <- coord_full[coord_full$Set != "Set average",]


### 3D Surface Area data ###

# Import file with all 3D areas
area_full <- read.csv(file.path(getwd(), "SuppInfo_5.csv"), 
                      header=TRUE, 
                      stringsAsFactors=FALSE, 
                      sep = ',')

# Create intra-observer dataframe
area_intra <- area_full[,-7]

# Create inter-observer dataframe
rowMeans <- rowMeans(area_full[,2:6])
area_inter <- data.frame(Name = area_full[,1], Obs1 = rowMeans, Obs2 = area_full[,7])


# Landmarks: Euclidean distance (pairs) -----------------------------------

### Intra-observer calculations ###

## Prepare data
# Separate datasets by landmark
coord_intra_LM1 <- coord_intra[,1:5]
names(coord_intra_LM1) <- c("Name", "Set", "x", "y", "z")
coord_intra_LM2 <- coord_intra[,c(1:2, 6:8)]
names(coord_intra_LM2) <- c("Name", "Set", "x", "y", "z")
coord_intra_LM3 <- coord_intra[,c(1:2, 9:11)]
names(coord_intra_LM3) <- c("Name", "Set", "x", "y", "z")


## Landmark 1
# Get all unique models
models <- unique(coord_intra$Name)

# Initialize a list to store results
intra_dist_LM1 <- list()

# Calculate the Euclidean distance between sets
for (model in models) {
  # Filter data for the current model
  model_data <- coord_intra_LM1 %>% filter(Name == model)
  
  # Get all unique sets and combinations
  sets <- unique(model_data$Set)
  set_combinations <- combn(sets, 2, simplify = FALSE)
  
  # Calculate distance for each combination
  distance <- lapply(set_combinations, function(sets) {
    set1_coords <- model_data %>% filter(Set == sets[1]) %>% select(x, y, z) %>% unlist()
    set2_coords <- model_data %>% filter(Set == sets[2]) %>% select(x, y, z) %>% unlist()
    
    dist <- calculate_euclidean_distance(set1_coords, set2_coords)
    
    # Return a named vector for easier conversion to data frame
    return(c(Model = model, Set1 = sets[1], Set2 = sets[2], Distance = round(dist, 2)))
  })
  
  # Convert results to data frame and store in the list
  intra_dist_LM1[[model]] <- as.data.frame(do.call(rbind, distance), stringsAsFactors = FALSE)
}

# Combine all results into a single data frame
intra_dist_LM1 <- do.call(rbind, intra_dist_LM1)
# Pivot to wide format
intra_dist_LM1 <- intra_dist_LM1 %>%
  unite("Set_Comparison", Set1, Set2, sep = " vs ") %>%
  pivot_wider(names_from = Set_Comparison, values_from = Distance)

# Calculate Mean and standard deviation
intra_dist_LM1 <- calculate_mean_sd(intra_dist_LM1)
intra_dist_LM1 <- append_overall_row(intra_dist_LM1)

print(intra_dist_LM1)


## Landmark 2
# Get all unique models
models <- unique(coord_intra$Name)

# Initialize a list to store results
intra_dist_LM2 <- list()

# Calculate the Euclidean distance between sets
for (model in models) {
  # Filter data for the current model
  model_data <- coord_intra_LM2 %>% filter(Name == model)
  
  # Get all unique sets and combinations
  sets <- unique(model_data$Set)
  set_combinations <- combn(sets, 2, simplify = FALSE)
  
  # Calculate distance for each combination
  distance <- lapply(set_combinations, function(sets) {
    set1_coords <- model_data %>% filter(Set == sets[1]) %>% select(x, y, z) %>% unlist()
    set2_coords <- model_data %>% filter(Set == sets[2]) %>% select(x, y, z) %>% unlist()
    
    dist <- calculate_euclidean_distance(set1_coords, set2_coords)
    
    # Return a named vector for easier conversion to data frame
    return(c(Model = model, Set1 = sets[1], Set2 = sets[2], Distance = round(dist, 2)))
  })
  
  # Convert results to data frame and store in the list
  intra_dist_LM2[[model]] <- as.data.frame(do.call(rbind, distance), stringsAsFactors = FALSE)
}

# Combine all results into a single data frame
intra_dist_LM2 <- do.call(rbind, intra_dist_LM2)
# Pivot to wide format
intra_dist_LM2 <- intra_dist_LM2 %>%
  unite("Set_Comparison", Set1, Set2, sep = " vs ") %>%
  pivot_wider(names_from = Set_Comparison, values_from = Distance)

# Calculate Mean and Standard Deviation
intra_dist_LM2 <- calculate_mean_sd(intra_dist_LM2)
intra_dist_LM2 <- append_overall_row(intra_dist_LM2)

print(intra_dist_LM2)


## Landmark 3
# Get all unique models
models <- unique(coord_intra$Name)

# Initialize a list to store results
intra_dist_LM3 <- list()

# Calculate the Euclidean distance between sets
for (model in models) {
  # Filter data for the current model
  model_data <- coord_intra_LM3 %>% filter(Name == model)
  
  # Get all unique sets and combinations
  sets <- unique(model_data$Set)
  set_combinations <- combn(sets, 2, simplify = FALSE)
  
  # Calculate distance for each combination
  distance <- lapply(set_combinations, function(sets) {
    set1_coords <- model_data %>% filter(Set == sets[1]) %>% select(x, y, z) %>% unlist()
    set2_coords <- model_data %>% filter(Set == sets[2]) %>% select(x, y, z) %>% unlist()
    
    dist <- calculate_euclidean_distance(set1_coords, set2_coords)
    
    # Return a named vector for easier conversion to data frame
    return(c(Model = model, Set1 = sets[1], Set2 = sets[2], Distance = round(dist, 2)))
  })
  
  # Convert results to data frame and store in the list
  intra_dist_LM3[[model]] <- as.data.frame(do.call(rbind, distance), stringsAsFactors = FALSE)
}

# Combine all results into a single data frame
intra_dist_LM3 <- do.call(rbind, intra_dist_LM3)
# Pivot to wide format
intra_dist_LM3 <- intra_dist_LM3 %>%
  unite("Set_Comparison", Set1, Set2, sep = " vs ") %>%
  pivot_wider(names_from = Set_Comparison, values_from = Distance)

# Calculate Mean and Standard Deviation
intra_dist_LM3 <- calculate_mean_sd(intra_dist_LM3)
intra_dist_LM3 <- append_overall_row(intra_dist_LM3)

print(intra_dist_LM3)


### Inter-observer calculations ###

## Prepare data
coord_inter_LM1 <- coord_full[,1:6]
names(coord_inter_LM1) <- c("Name", "Observer", "Set", "x", "y", "z")
coord_inter_LM1_clean <- coord_inter_LM1[coord_inter_LM1$Name != "Model 20",]
coord_inter_LM2 <- coord_full[,c(1:3, 7:9)]
names(coord_inter_LM2) <- c("Name", "Observer", "Set", "x", "y", "z")
coord_inter_LM3 <- coord_full[,c(1:3, 10:12)]
names(coord_inter_LM3) <- c("Name", "Observer", "Set", "x", "y", "z")


## Landmark 1
# Get all unique models
models <- unique(coord_inter_LM1_clean$Name)

# Initialize a list to store results
inter_dist_LM1 <- list()

# Calculate Euclidean distance between dataset of Observer 2 and each set of Observer 1
for (model in models) {
  # Filter data for the current model
  model_data <- coord_inter_LM1_clean %>% filter(Name == model)
  
  # Get Observer 1 and Observer 2 data
  observer1_data <- model_data %>% filter(Observer == 1)
  observer2_data <- model_data %>% filter(Observer == 2)
  
  # Extract the unique sets for Observer 1
  sets <- unique(observer1_data$Set)
  
  # Calculate distance between Observer 2 and each set of Observer 1
  distance <- sapply(sets, function(set1) {
    dist <- sqrt(
      (observer1_data$x[observer1_data$Set == set1] - observer2_data$x)^2 +
        (observer1_data$y[observer1_data$Set == set1] - observer2_data$y)^2 +
        (observer1_data$z[observer1_data$Set == set1] - observer2_data$z)^2)
    
    return(c(Model = model, Set1 = set1, Set2 = "Observer 2", Distance = round(dist, 2)))
  })
  
  # Convert results to data frame and add to list
  inter_dist_LM1[[model]] <- as.data.frame(t(distance), stringsAsFactors = FALSE)
}

# Combine all results into a single data frame
inter_dist_LM1 <- do.call(rbind, inter_dist_LM1)

# Pivot to wide format
inter_dist_LM1 <- inter_dist_LM1 %>%
  unite("Set_Comparison", Set1, Set2, sep = " vs ") %>%
  pivot_wider(names_from = Set_Comparison, values_from = Distance)

# Calculate Mean and Standard Deviation
inter_dist_LM1 <- calculate_mean_sd(inter_dist_LM1)
inter_dist_LM1 <- append_overall_row(inter_dist_LM1)

print(inter_dist_LM1)


## Landmark 2
# Get all unique models
models <- unique(coord_inter_LM2$Name)
# Initialize a list to store results
inter_dist_LM2 <- list()

# Calculate Euclidean distance between dataset of Observer 2 and each set of Observer 1
for (model in models) {
  # Filter data for the current model
  model_data <- coord_inter_LM2 %>% filter(Name == model)
  
  # Get Observer 1 and Observer 2 data
  observer1_data <- model_data %>% filter(Observer == 1)
  observer2_data <- model_data %>% filter(Observer == 2)
  
  # Extract the unique sets for Observer 1
  sets <- unique(observer1_data$Set)
  
  # Calculate distance between Observer 2 and each set of Observer 1
  distance <- sapply(sets, function(set1) {
    dist <- sqrt(
      (observer1_data$x[observer1_data$Set == set1] - observer2_data$x)^2 +
        (observer1_data$y[observer1_data$Set == set1] - observer2_data$y)^2 +
        (observer1_data$z[observer1_data$Set == set1] - observer2_data$z)^2)
    
    return(c(Model = model, Set1 = set1, Set2 = "Observer 2", Distance = round(dist, 2)))
  })
  
  # Convert results to data frame and add to list
  inter_dist_LM2[[model]] <- as.data.frame(t(distance), stringsAsFactors = FALSE)
}

# Combine all results into a single data frame
inter_dist_LM2 <- do.call(rbind, inter_dist_LM2)
# Pivot to wide format
inter_dist_LM2 <- inter_dist_LM2 %>%
  unite("Set_Comparison", Set1, Set2, sep = " vs ") %>%
  pivot_wider(names_from = Set_Comparison, values_from = Distance)

# Calculate Mean and Standard Deviation
inter_dist_LM2 <- calculate_mean_sd(inter_dist_LM2)
inter_dist_LM2 <- append_overall_row(inter_dist_LM2)

print(inter_dist_LM2)


## Landmark 3
# Get all unique models
models <- unique(coord_inter_LM3$Name)
# Initialize a list to store results
inter_dist_LM3 <- list()

# Calculate Euclidean distance between dataset of Observer 2 and each set of Observer 1
for (model in models) {
  # Filter data for the current model
  model_data <- coord_inter_LM3 %>% filter(Name == model)
  
  # Get Observer 1 and Observer 2 data
  observer1_data <- model_data %>% filter(Observer == 1)
  observer2_data <- model_data %>% filter(Observer == 2)
  
  # Extract the unique sets for Observer 1
  sets <- unique(observer1_data$Set)
  
  # Calculate distance between Observer 2 and each set of Observer 1
  distance <- sapply(sets, function(set1) {
    dist <- sqrt(
      (observer1_data$x[observer1_data$Set == set1] - observer2_data$x)^2 +
        (observer1_data$y[observer1_data$Set == set1] - observer2_data$y)^2 +
        (observer1_data$z[observer1_data$Set == set1] - observer2_data$z)^2)
    
    return(c(Model = model, Set1 = set1, Set2 = "Observer 2", Distance = round(dist, 2)))
  })
  
  # Convert results to data frame and add to list
  inter_dist_LM3[[model]] <- as.data.frame(t(distance), stringsAsFactors = FALSE)
}

# Combine all results into a single data frame
inter_dist_LM3 <- do.call(rbind, inter_dist_LM3)
# Pivot to wide format
inter_dist_LM3 <- inter_dist_LM3 %>%
  unite("Set_Comparison", Set1, Set2, sep = " vs ") %>%
  pivot_wider(names_from = Set_Comparison, values_from = Distance)

# Calculate Mean and Standard Deviation
inter_dist_LM3 <- calculate_mean_sd(inter_dist_LM3)
inter_dist_LM3 <- append_overall_row(inter_dist_LM3)

print(inter_dist_LM3)


### Present data as boxplot ###

## Reshape and combine all datasets
# Function to reshape data
reshape_data <- function(data, observer_type, landmark_id) {
  data %>%
    pivot_longer(cols = -Model, names_to = "Comparison", values_to = "Value") %>%
    mutate(Observer = observer_type, Landmark = landmark_id)
}

# Reshape all intra observer datasets
intra_LM1_long <- reshape_data(intra_dist_LM1[,1:11], "Intra-observer", "1")
intra_LM2_long <- reshape_data(intra_dist_LM2[,1:11], "Intra-observer", "2")
intra_LM3_long <- reshape_data(intra_dist_LM3[,1:11], "Intra-observer", "3")

# Reshape all inter observer datasets
inter_LM1_long <- reshape_data(inter_dist_LM1[,1:6], "Inter-observer", "1")
inter_LM2_long <- reshape_data(inter_dist_LM2[,1:6], "Inter-observer", "2")
inter_LM3_long <- reshape_data(inter_dist_LM3[,1:6], "Inter-observer", "3")

# Combine all reshaped data
combined_data <- bind_rows(
  intra_LM1_long,
  intra_LM2_long,
  intra_LM3_long,
  inter_LM1_long,
  inter_LM2_long,
  inter_LM3_long
)

# Switch the order of the Observer levels to have "Intra-observer" first
combined_data$Observer <- factor(combined_data$Observer, 
                                 levels = c("Intra-observer", "Inter-observer"))

# Create the boxplot
dist_plot <- ggplot(combined_data, aes(x = Landmark, y = Value, fill = Observer)) +
  geom_boxplot(outlier.alpha = 0.8,
               outlier.size = 1, 
               outlier.shape = 16) + 
  labs(x = "Landmark",
       y = "Euclidean distance (mm)",
       fill = "Comparison") +
  theme_minimal() +
  scale_fill_manual(values = c("Intra-observer" = "#FBE49A", "Inter-observer" = "#F28400")) +  
  theme(text = element_text(family = "sans"),
        axis.text.x = element_text(hjust = 1, vjust = 1, size = 10),  
        axis.title.x = element_text(margin = margin(t = 10), size = 13),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(r = 10), size = 13),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 13),
        legend.key.size = unit(1.5, "lines"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = c(0.98, 0.97),  
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "#FFFFFF", color = "#D3D3D3"))

dist_plot


# Landmarks: Distance-based ICC -------------------------------------------

### Intra-observer calculations ###

# Prepare data frame for results
coord_intra_dICC <- data.frame(
  LM = character(),
  ICC = numeric(),
  SE = numeric(),
  Lower_CI = numeric(),
  Upper_CI = numeric(),
  Iterations = integer(),  # Add an Iterations column
  stringsAsFactors = FALSE
)

# Specify the models
strata <- factor(coord_intra$Name)

# Run the analysis for the different landmarks
coord_intra_dICC <- calculate_icc(coord_intra, "Intra_LM1", "LM1_x", "LM1_y", "LM1_z", strata, coord_intra_dICC)
coord_intra_dICC <- calculate_icc(coord_intra, "Intra_LM2", "LM2_x", "LM2_y", "LM2_z", strata, coord_intra_dICC)
coord_intra_dICC <- calculate_icc(coord_intra, "Intra_LM3", "LM3_x", "LM3_y", "LM3_z", strata, coord_intra_dICC)

print(coord_intra_dICC)


### Inter-observer calculations ###

## Prepare data frame for results
coord_inter_dICC <- data.frame(
  LM = character(),
  ICC = numeric(),
  SE = numeric(),
  Lower_CI = numeric(),
  Upper_CI = numeric(),
  Iterations = integer(),  # Add an Iterations column
  stringsAsFactors = FALSE
)

## Landmark 1
# Specify the models and clean up dataset (removal of Model 20)
coord_inter_clean <- coord_inter[coord_inter$Name != "Model 20",] 
strata <- factor(coord_inter_clean$Name)

# Run the analysis
coord_inter_dICC <- calculate_icc(coord_inter_clean, "Inter_LM1", "LM1_x", 
                                  "LM1_y", "LM1_z", strata, coord_inter_dICC)
coord_full_clean <- coord_full[coord_full$Name != "Model 20",]
temp <- coord_full_clean[coord_full_clean$Set == "Set A",]
coord_inter_dICC <- calculate_icc(temp, "Inter_LM1_1Avs2A", "LM1_x", "LM1_y", 
                                  "LM1_z", strata, coord_inter_dICC)
temp <- coord_full_clean[coord_full_clean$Set == "Set B" | coord_full_clean$Observer == 2,]
coord_inter_dICC <- calculate_icc(temp, "Inter_LM1_1Bvs2A", "LM1_x", "LM1_y", 
                                  "LM1_z", strata, coord_inter_dICC)
temp <- coord_full_clean[coord_full_clean$Set == "Set C" | coord_full_clean$Observer == 2,]
coord_inter_dICC <- calculate_icc(temp, "Inter_LM1_1Cvs2A", "LM1_x", "LM1_y", 
                                  "LM1_z", strata, coord_inter_dICC)
temp <- coord_full_clean[coord_full_clean$Set == "Set D" | coord_full_clean$Observer == 2,]
coord_inter_dICC <- calculate_icc(temp, "Inter_LM1_1Dvs2A", "LM1_x", "LM1_y", 
                                  "LM1_z", strata, coord_inter_dICC)
temp <- coord_full_clean[coord_full_clean$Set == "Set E" | coord_full_clean$Observer == 2,]
coord_inter_dICC <- calculate_icc(temp, "Inter_LM1_1Evs2A", "LM1_x", "LM1_y", 
                                  "LM1_z", strata, coord_inter_dICC)


## Landmarks 2 and 3
strata <- factor(coord_inter$Name)
coord_inter_dICC <- calculate_icc(coord_inter, "Inter_LM2", "LM2_x", "LM2_y", 
                                  "LM2_z", strata, coord_inter_dICC)
coord_inter_dICC <- calculate_icc(coord_inter, "Inter_LM3", "LM3_x", "LM3_y", 
                                  "LM3_z", strata, coord_inter_dICC)

temp <- coord_full[coord_full$Set == "Set A",]
coord_inter_dICC <- calculate_icc(temp, "Inter_LM2_1Avs2A", "LM2_x", "LM2_y", 
                                  "LM2_z", strata, coord_inter_dICC)
coord_inter_dICC <- calculate_icc(temp, "Inter_LM3_1Avs2A", "LM3_x", "LM3_y", 
                                  "LM3_z", strata, coord_inter_dICC)
temp <- coord_full[coord_full$Set == "Set B" | coord_full$Observer == 2,]
coord_inter_dICC <- calculate_icc(temp, "Inter_LM2_1Bvs2A", "LM2_x", "LM2_y", 
                                  "LM2_z", strata, coord_inter_dICC)
coord_inter_dICC <- calculate_icc(temp, "Inter_LM3_1Bvs2A", "LM3_x", "LM3_y", 
                                  "LM3_z", strata, coord_inter_dICC)
temp <- coord_full[coord_full$Set == "Set C" | coord_full$Observer == 2,]
coord_inter_dICC <- calculate_icc(temp, "Inter_LM2_1Cvs2A", "LM2_x", "LM2_y", 
                                  "LM2_z", strata, coord_inter_dICC)
coord_inter_dICC <- calculate_icc(temp, "Inter_LM3_1Cvs2A", "LM3_x", "LM3_y", 
                                  "LM3_z", strata, coord_inter_dICC)
temp <- coord_full[coord_full$Set == "Set D" | coord_full$Observer == 2,]
coord_inter_dICC <- calculate_icc(temp, "Inter_LM2_1Dvs2A", "LM2_x", "LM2_y", 
                                  "LM2_z", strata, coord_inter_dICC)
coord_inter_dICC <- calculate_icc(temp, "Inter_LM3_1Dvs2A", "LM3_x", "LM3_y", 
                                  "LM3_z", strata, coord_inter_dICC)
temp <- coord_full[coord_full$Set == "Set E" | coord_full$Observer == 2,]
coord_inter_dICC <- calculate_icc(temp, "Inter_LM2_1Evs2A", "LM2_x", "LM2_y", 
                                  "LM2_z", strata, coord_inter_dICC)
coord_inter_dICC <- calculate_icc(temp, "Inter_LM3_1Evs2A", "LM3_x", "LM3_y", 
                                  "LM3_z", strata, coord_inter_dICC)

print(coord_inter_dICC)

# Landmark data: Export statistical results -------------------------------

### Automatically export the boxplot as a PDF file ###
ggsave("boxplot_eucl_dist.pdf", plot = dist_plot, width = 200, height = 130, 
       units = "mm", dpi = 1000)


### Export distance results as excel file ###
# Create a new workbook
wb <- createWorkbook()

# Write data frames to worksheets
addWorksheet(wb, "Intra-observer LM1")
writeData(wb, "Intra-observer LM1", intra_dist_LM1)

addWorksheet(wb, "Intra-observer LM2")
writeData(wb, "Intra-observer LM2", intra_dist_LM2)

addWorksheet(wb, "Intra-observer LM3")
writeData(wb, "Intra-observer LM3", intra_dist_LM3)

addWorksheet(wb, "Inter-observer LM1")
writeData(wb, "Inter-observer LM1", inter_dist_LM1)

addWorksheet(wb, "Inter-observer LM2")
writeData(wb, "Inter-observer LM2", inter_dist_LM2)

addWorksheet(wb, "Inter-observer LM3")
writeData(wb, "Inter-observer LM3", inter_dist_LM3)

# Save the workbook
saveWorkbook(wb, file = file.path(getwd(), "LandmarkCoordinates_Distance(Pairs).xlsx"), 
             overwrite = TRUE)


### Export dICC results as excel file ###
# Create a new workbook
wb <- createWorkbook()

# Write data frames to worksheets
addWorksheet(wb, "Intra-observer dICC")
writeData(wb, "Intra-observer dICC", coord_intra_dICC)

addWorksheet(wb, "Inter-observer dICC")
writeData(wb, "Inter-observer dICC", coord_inter_dICC)

# Save the workbook
saveWorkbook(wb, file = file.path(getwd(), "LandmarkCoordinates_dICC.xlsx"), 
             overwrite = TRUE)


# 3D Area: descriptive statistics -----------------------------------------

### Intra-observer calculations ###

# Calculate mean and standard deviation for each model
rowmean <- rowMeans(area_intra[,-1], na.rm = TRUE)
row_stdev <- apply(area_intra[,-1], 1, sd, na.rm = TRUE)

# Calculate mean and standard deviation for each set
colmean <- colMeans(area_intra[,-1], na.rm = TRUE)
col_stdev <- sapply(area_intra[,-1], sd, na.rm = TRUE)

# Calculate overall mean and standard deviation
mean_overall <- mean(unlist(area_intra[,-1]), na.rm = TRUE)
SD_overall <- sd(unlist(area_intra[,-1]), na.rm = TRUE)

# Construct data frame
area_intra_descr <- area_intra
area_intra_descr$Mean <- rowmean
area_intra_descr$SD <- row_stdev
area_intra_descr$CV <- (row_stdev/rowmean)*100
area_intra_descr <- rbind(area_intra_descr, c("Mean", colmean, mean_overall, NA, NA))
area_intra_descr <- rbind(area_intra_descr, c("SD", col_stdev, NA, SD_overall, NA))
area_intra_descr <- rbind(area_intra_descr, c("CV", (col_stdev/colmean)*100, NA, NA, (SD_overall/mean_overall)*100))


### Inter-observer calculations ###

# Removal of model 20 in Observer 1
area_inter_clean <- na.omit(area_inter)

# Calculate mean and Standard Deviation for each set
colmean <- colMeans(area_inter_clean[,-1], na.rm = TRUE)
col_stdev <- sapply(area_inter_clean[,-1], sd, na.rm = TRUE)

# Construct data frame
area_inter_descr <- rbind(area_inter_clean, c("Mean", colmean))
area_inter_descr <- rbind(area_inter_descr, c("SD", col_stdev))
area_inter_descr <- rbind(area_inter_descr, c("CV", (col_stdev/colmean)*100))


# 3D Area: percentage error -----------------------------------------------

### Intra-observer calculations ###

# Calculate percentage error between sets
area_intra_perc_AvsB <- percent_error(area_intra$Obs1_SetA, area_intra$Obs1_SetB)
area_intra_perc_AvsC <- percent_error(area_intra$Obs1_SetA, area_intra$Obs1_SetC)
area_intra_perc_AvsD <- percent_error(area_intra$Obs1_SetA, area_intra$Obs1_SetD)
area_intra_perc_AvsE <- percent_error(area_intra$Obs1_SetA, area_intra$Obs1_SetE)
area_intra_perc_BvsC <- percent_error(area_intra$Obs1_SetB, area_intra$Obs1_SetC)
area_intra_perc_BvsD <- percent_error(area_intra$Obs1_SetB, area_intra$Obs1_SetD)
area_intra_perc_BvsE <- percent_error(area_intra$Obs1_SetB, area_intra$Obs1_SetE)
area_intra_perc_CvsD <- percent_error(area_intra$Obs1_SetC, area_intra$Obs1_SetD)
area_intra_perc_CvsE <- percent_error(area_intra$Obs1_SetC, area_intra$Obs1_SetE)
area_intra_perc_DvsE <- percent_error(area_intra$Obs1_SetD, area_intra$Obs1_SetE)

# Construct data frame
area_intra_perc <- data.frame(Name = area_intra$Name, AvsB = area_intra_perc_AvsB, 
                              AvsC= area_intra_perc_AvsC, AvsD = area_intra_perc_AvsD, 
                              AvsE = area_intra_perc_AvsE, BvsC = area_intra_perc_BvsC, 
                              BvsD = area_intra_perc_BvsD, BvsE = area_intra_perc_BvsE, 
                              CvsD = area_intra_perc_CvsD, CvsE = area_intra_perc_CvsE, 
                              DvsE = area_intra_perc_DvsE)

# Calculate mean and standard deviation in percentage error for each model
rowmean <- rowMeans(area_intra_perc[,-1], na.rm = TRUE)
row_stdev <- apply(area_intra_perc[,-1], 1, sd, na.rm = TRUE)

# Calculate mean and standard deviation in percentage error for each set
colmean <- colMeans(area_intra_perc[,-1], na.rm = TRUE)
col_stdev <- sapply(area_intra_perc[,-1], sd, na.rm = TRUE)

# Calculate overall mean and standard deviation in percentage error
mean_overall <- mean(unlist(area_intra_perc[,-1]), na.rm = TRUE)
SD_overall <- sd(unlist(area_intra_perc[,-1]), na.rm = TRUE)

# Add descriptive statistics to dataframe
area_intra_perc$Average <- rowmean
area_intra_perc$SD <- row_stdev
area_intra_perc <- rbind(area_intra_perc, c("Average Percentage Error", colmean, mean_overall, NA))
area_intra_perc <- rbind(area_intra_perc, c("SD", col_stdev, NA, SD_overall))

# Round numeric values for export
area_intra_perc <- area_intra_perc %>%
  mutate(across(-Name, as.numeric)) %>% 
  mutate(across(-Name, round, digits = 2))


### Inter-observer calculations ###

## Percentage calculations for all sets of Observer 1 and Observer 2
# Clean up dataset (removal of Model 20)
area_full_clean <- area_full[area_full$Name != "Model 20",]

# Calculate percentage error between sets
area_inter_perc_1Avs2A <- percent_error(area_full_clean$Obs1_SetA, area_full_clean$Obs2_SetA)
area_inter_perc_1Bvs2A <- percent_error(area_full_clean$Obs1_SetB, area_full_clean$Obs2_SetA)
area_inter_perc_1Cvs2A <- percent_error(area_full_clean$Obs1_SetC, area_full_clean$Obs2_SetA)
area_inter_perc_1Dvs2A <- percent_error(area_full_clean$Obs1_SetD, area_full_clean$Obs2_SetA)
area_inter_perc_1Evs2A <- percent_error(area_full_clean$Obs1_SetE, area_full_clean$Obs2_SetA)

# Construct data frame
area_inter_perc <- data.frame(Name = area_full_clean$Name, 
                              Avs2A = area_inter_perc_1Avs2A, 
                              Bvs2A= area_inter_perc_1Bvs2A,
                              Cvs2A = area_inter_perc_1Cvs2A, 
                              Dvs2A = area_inter_perc_1Dvs2A,
                              Evs2A = area_inter_perc_1Evs2A)

# Calculate mean and standard deviation for each model
rowmean <- rowMeans(area_inter_perc[,-1], na.rm = TRUE)
row_stdev <- apply(area_inter_perc[,-1], 1, sd, na.rm = TRUE)

# Calculate mean and standard deviation for each set
colmean <- colMeans(area_inter_perc[,-1], na.rm = TRUE)
col_stdev <- sapply(area_inter_perc[,-1], sd, na.rm = TRUE)

# Calculate overall mean and standard deviation
mean_overall <- mean(unlist(area_inter_perc[,-1]), na.rm = TRUE)
SD_overall <- sd(unlist(area_inter_perc[,-1]), na.rm = TRUE)

# Add descriptive statistics to dataframe
area_inter_perc$Average <- rowmean
area_inter_perc$SD <- row_stdev
area_inter_perc <- rbind(area_inter_perc, c("Average Percentage Error", colmean, mean_overall, NA))
area_inter_perc <- rbind(area_inter_perc, c("SD", col_stdev, NA, SD_overall))

# Round numeric values for export
area_inter_perc <- area_inter_perc %>%
  mutate(across(-Name, as.numeric)) %>% 
  mutate(across(-Name, round, digits = 2))


## Percentage calculations with average of Observer 1 and Observer 2
# Clean up dataset (removal of Model 20)
area_inter_perc_av <- area_inter[area_inter$Name != "Model 20",]

# Calculate percent errors for each observation
area_inter_perc_av$Percent_Error <- percent_error(area_inter_perc_av$Obs1, area_inter_perc_av$Obs2)

# Calculate mean and standard deviation for the percentage column
colmean <- mean(area_inter_perc_av$Percent_Error, na.rm = TRUE)
col_stdev <- sd(area_inter_perc_av$Percent_Error, na.rm = TRUE)
area_inter_perc_av <- rbind(area_inter_perc_av, c("Average Percentage Error", NA, NA, colmean))
area_inter_perc_av <- rbind(area_inter_perc_av, c("SD Percentage Error", NA, NA, col_stdev))

# Round numeric values for export
area_inter_perc_av <- area_inter_perc_av %>%
  mutate(across(-Name, as.numeric)) %>% 
  mutate(across(-Name, round, digits = 2))


# 3D Area: correlation coefficients ---------------------------------------

### Test for normality ###

# Data is normally distributed
area_intra_A <- shapiro.test(area_intra[,2])
area_intra_B <- shapiro.test(area_intra[,3])
area_intra_C <- shapiro.test(area_intra[,4])
area_intra_D <- shapiro.test(area_intra[,5])
area_intra_E <- shapiro.test(area_intra[,6])
area_inter_obs1 <- shapiro.test(area_inter[,2])
area_inter_obs2 <- shapiro.test(area_inter[,3])

# Export shapiro-wilk's results
area_normality <- data.frame(Name = c("Set 1A", "Set 1B", "Set 1C", "Set 1D", 
                                      "Set 1E", "Av. obs 1", "Set 2A"),
                             W = c(area_intra_A$statistic, area_intra_B$statistic, 
                                   area_intra_C$statistic, area_intra_D$statistic, 
                                   area_intra_E$statistic, area_inter_obs1$statistic,
                                   area_inter_obs2$statistic),
                             p.value = c(area_intra_A$p.value, area_intra_B$p.value, 
                                         area_intra_C$p.value, area_intra_D$p.value, 
                                         area_intra_E$p.value, area_inter_obs1$p.value, 
                                         area_inter_obs2$p.value))


### Intra-observer calculations ###

# Calculate ICC
area_intra_icc <- icc(area_intra[,-1], 
                      model = "twoway", 
                      type = "agreement", 
                      unit = "single")

# Calculate Lin's CCC
area_intra_ccc_AvsB <- CCC(area_intra[,2], area_intra[,3])
area_intra_ccc_AvsC <- CCC(area_intra[,2], area_intra[,4])
area_intra_ccc_AvsD <- CCC(area_intra[,2], area_intra[,5])
area_intra_ccc_AvsE <- CCC(area_intra[,2], area_intra[,6])
area_intra_ccc_BvsC <- CCC(area_intra[,3], area_intra[,4])
area_intra_ccc_BvsD <- CCC(area_intra[,3], area_intra[,5])
area_intra_ccc_BvsE <- CCC(area_intra[,3], area_intra[,6])
area_intra_ccc_CvsD <- CCC(area_intra[,4], area_intra[,5])
area_intra_ccc_CvsE <- CCC(area_intra[,4], area_intra[,6])
area_intra_ccc_DvsE <- CCC(area_intra[,5], area_intra[,6])


### Inter-observer calculations ###

# Calculate ICC
area_intra_clean <- area_intra[area_intra$Name != "Model 20",]
area_clean <- merge(area_intra_clean, area_inter_clean, by="Name")
area_inter_icc <- icc(area_inter_clean[,-1], model = "twoway", 
                      type = "agreement", unit = "single")
area_inter_icc_1Avs2A <- icc(area_clean[,c(2,8)], model = "twoway", 
                             type = "agreement", unit = "single")
area_inter_icc_1Bvs2A <- icc(area_clean[,c(3,8)], model = "twoway", 
                             type = "agreement", unit = "single")
area_inter_icc_1Cvs2A <- icc(area_clean[,c(4,8)], model = "twoway", 
                             type = "agreement", unit = "single")
area_inter_icc_1Dvs2A <- icc(area_clean[,c(5,8)], model = "twoway", 
                             type = "agreement", unit = "single")
area_inter_icc_1Evs2A <- icc(area_clean[,c(6,8)], model = "twoway", 
                             type = "agreement", unit = "single")

# Calculate Lin's CCC
area_inter_ccc <- CCC(area_inter_clean[,2], area_inter_clean[,3])
area_inter_ccc_1Avs2A <- CCC(area_intra_clean[,2], area_inter_clean[,3])
area_inter_ccc_1Bvs2A <- CCC(area_intra_clean[,3], area_inter_clean[,3])
area_inter_ccc_1Cvs2A <- CCC(area_intra_clean[,4], area_inter_clean[,3])
area_inter_ccc_1Dvs2A <- CCC(area_intra_clean[,5], area_inter_clean[,3])
area_inter_ccc_1Evs2A <- CCC(area_intra_clean[,6], area_inter_clean[,3])


### Prepare results for export ###

## ICC
# Function to extract ICC components and return as a row
extract_icc_results <- function(name, icc_result) {
  data.frame(
    Name = name,
    Sample = icc_result$subjects,
    Observers = icc_result$raters,
    ICC = icc_result$value,
    LB = icc_result$lbound,
    UB = icc_result$ubound,
    F.value = icc_result$Fvalue,
    F.df1 = icc_result$df1,
    F.df2 = icc_result$df2,
    Sig = icc_result$p.value,
    stringsAsFactors = FALSE
  )
}

# List of comparison names and corresponding CCC objects
comparisons <- list(
  "Intra-observer" = area_intra_icc,
  "Inter-observer" = area_inter_icc,
  "1A vs. 2A" = area_inter_icc_1Avs2A,
  "1B vs. 2A" = area_inter_icc_1Bvs2A,
  "1C vs. 2A" = area_inter_icc_1Cvs2A,
  "1D vs. 2A" = area_inter_icc_1Dvs2A,
  "1E vs. 2A" = area_inter_icc_1Evs2A
)

# Initialize an empty data frame to store results
area_icc_result <- data.frame()

# Loop through the comparisons and extract ICC results
for (name in names(comparisons)) {
  area_icc_result <- rbind(area_icc_result, 
                           extract_icc_results(name, comparisons[[name]]))
}

area_icc_result <- area_icc_result %>%
  mutate(across(where(is.numeric), round, digits = 3))

print(area_icc_result)


## Lin's CCC
# Function to extract CCC components and return as a row
extract_ccc_results <- function(name, ccc_result) {
  data.frame(
    Name = name,
    CCC = ccc_result$rho.c$est,
    LB = ccc_result$rho.c$lwr.ci,
    UB = ccc_result$rho.c$upr.ci,
    Bias = ccc_result$C.b,
    Scale = ccc_result$s.shift,
    Location = ccc_result$l.shift,
    stringsAsFactors = FALSE
  )
}

# List of comparison names and corresponding CCC objects
comparisons <- list(
  "1A vs. 1B" = area_intra_ccc_AvsB,
  "1A vs. 1C" = area_intra_ccc_AvsC,
  "1A vs. 1D" = area_intra_ccc_AvsD,
  "1A vs. 1E" = area_intra_ccc_AvsE,
  "1B vs. 1C" = area_intra_ccc_BvsC,
  "1B vs. 1D" = area_intra_ccc_BvsD,
  "1B vs. 1E" = area_intra_ccc_BvsE,
  "1C vs. 1D" = area_intra_ccc_CvsD,
  "1C vs. 1E" = area_intra_ccc_CvsE,
  "1D vs. 1E" = area_intra_ccc_DvsE,
  "Inter-observer" = area_inter_ccc,
  "1A vs. 2A" = area_inter_ccc_1Avs2A,
  "1B vs. 2A" = area_inter_ccc_1Bvs2A,
  "1C vs. 2A" = area_inter_ccc_1Cvs2A,
  "1D vs. 2A" = area_inter_ccc_1Dvs2A,
  "1E vs. 2A" = area_inter_ccc_1Evs2A
)

# Initialize an empty data frame to store results
area_ccc_result <- data.frame()

# Loop through the comparisons and extract CCC results
for (name in names(comparisons)) {
  area_ccc_result <- rbind(area_ccc_result, 
                           extract_ccc_results(name, comparisons[[name]]))
}

area_ccc_result <- area_ccc_result %>%
  mutate(across(where(is.numeric), round, digits = 3))

print(area_ccc_result)


### Results of Lin's CCC presented as a correlation matrix ###

## Create correlation matrix
# Filter out the 'Inter-observer' row
area_ccc_result_filtered <- area_ccc_result %>%
  filter(Name != "Inter-observer")

# Split the 'Name' column into 'Variable1' and 'Variable2'
area_ccc_split <- area_ccc_result_filtered %>%
  separate(Name, into = c("Variable1", "Variable2"), sep = " vs. ")

# Create an empty correlation matrix
variables <- unique(c(area_ccc_split$Variable1, area_ccc_split$Variable2))
ccc_matrix <- matrix(NA, nrow = length(variables), ncol = length(variables))

# Assign row and column names
rownames(ccc_matrix) <- colnames(ccc_matrix) <- variables

# Fill the correlation matrix
for (i in 1:nrow(area_ccc_split)) {
  var1 <- area_ccc_split$Variable1[i]
  var2 <- area_ccc_split$Variable2[i]
  ccc_value <- area_ccc_split$CCC[i]
  
  ccc_matrix[var1, var2] <- ccc_value
  ccc_matrix[var2, var1] <- ccc_value
}

# Fill diagonal with 1s
diag(ccc_matrix) <- 1

print(ccc_matrix)


## Present matrix as heatmap
# Define custom color palette
coul <- colorRampPalette(brewer.pal(9, "YlOrRd")[2:8])(256)

# Convert the matrix to long format
ccc_matrix_lhalf <- ccc_matrix
ccc_matrix_lhalf[lower.tri(ccc_matrix_lhalf)] <- NA
ccc_long <- melt(ccc_matrix_lhalf, varnames = c("Var1", "Var2"), value.name = "Value")

# Assign labels to the factors
ccc_long$Var1 <- factor(ccc_long$Var1, levels = unique(ccc_long$Var1))
ccc_long$Var2 <- factor(ccc_long$Var2, levels = unique(ccc_long$Var2))
ccc_long$Var2 <- factor(ccc_long$Var2, levels = rev(unique(ccc_long$Var2)))

# Create heatmap with custom colors and set NA values to white
heatmap_plot <- ggplot(data = ccc_long, aes(x = Var1, y = Var2, fill = Value)) +
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(Value), "", format(round(Value, 3), nsmall = 3))), 
            color = "black", size = 4) +
  scale_fill_gradientn(colors = coul, 
                       values = scales::rescale(c(0, 0.4, 0.650, 0.9, 1)),
                       guide = "colorbar", 
                       name = "Coefficient",
                       na.value = "white") + 
  labs(title = "Pairwise calculation of Lin's CCC") +
  theme_minimal() +
  theme(text = element_text(family = "sans"),
        axis.text.x = element_text(hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18))

heatmap_plot


### Concordance correlation plots ###

## Intra-observer plots
# Define the combinations and their corresponding labels
combinations <- c("AvsB", "AvsC", "AvsD", "AvsE", "BvsC", "BvsD", "BvsE", 
                  "CvsD", "CvsE", "DvsE")

# Loop through each combination
for (comb in combinations) {
  ccc_variable <- get(paste0("area_intra_ccc_", comb))
  
  # Extract the sets from the combination string
  sets <- unlist(strsplit(comb, "vs"))
  
  # Ensure the area_intra variable is indexed correctly by column names
  var1 <- area_intra[[paste0("Obs1_Set", sets[1])]]
  var2 <- area_intra[[paste0("Obs1_Set", sets[2])]]
  
  # Create the label for the CCC value and CI
  lab <- paste("CCC: ", round(ccc_variable$rho.c[, 1], digits = 2), ", 95% CI [", 
               round(ccc_variable$rho.c[, 2], digits = 2), " - ", 
               round(ccc_variable$rho.c[, 3], digits = 2), "]", sep = "")
  
  # Create the PNG file name based on the combination
  png_filename <- paste0("concordance_correlation_plot_Intra_", comb, ".png")
  png(png_filename, width = 200, height = 200, units = "mm", res = 1000)
  
  # Plotting
  plot(var1, var2, 
       xlab = paste("Set", sets[1]), 
       ylab = paste("Set", sets[2]), 
       pch = 16,
       main = paste("Intra-observer comparison: Set", sets[1], "and Set", sets[2]),
       cex.main = 1.6, cex.lab = 1.3, cex.axis = 0.8)
  
  # Add line of perfect concordance
  abline(a = 0, b = 1, lty = 2, col = "red", lwd = 2)
  
  # Calculate RMA slope and intercept
  sd_x <- sd(var1)
  sd_y <- sd(var2)
  rma_slope <- sd_y / sd_x
  rma_intercept <- mean(var2) - rma_slope * mean(var1)
  
  # Add the RMA line
  abline(a = rma_intercept, b = rma_slope, lty = 1, lwd = 2)
  
  # Add legend
  legend("topleft", legend = c("Line of perfect concordance", "Reduced major axis"), 
         lty = c(2, 1), col = c("red", "black"), lwd = c(2, 2), bty = "n", cex = 1)
  
  # Position text in the lower right corner
  x_pos <- max(var1) * 0.8
  y_pos <- min(var2)
  
  # Add text
  text(x = x_pos, y = y_pos, labels = lab, pos = 4, cex = 0.9, col = "black")

  dev.off()
}


## Inter-observer plots
# Define the combinations and their corresponding labels
area_inter_ccc_1vs2A <- area_inter_ccc
names(area_clean) <- c("Name", "Set1A", "Set1B", "Set1C", "Set1D", "Set1E", "Set1", "Set2A")
combinations <- c("1Avs2A", "1Bvs2A", "1Cvs2A", "1Dvs2A", "1Evs2A", "1vs2A")

# Loop through each combination
for (comb in combinations) {
  # Dynamically create the CCC variable name
  ccc_variable <- get(paste0("area_inter_ccc_", comb))
  
  # Extract the sets from the combination string
  sets <- unlist(strsplit(comb, "vs"))
  
  # Ensure the area_intra variable is indexed correctly by column names
  var1 <- area_clean[[paste0("Set", sets[1])]]
  var2 <- area_clean[[paste0("Set", sets[2])]]
  
  # Create the label for the CCC value and CI
  lab <- paste("CCC: ", round(ccc_variable$rho.c[, 1], digits = 2), ", 95% CI [", 
               round(ccc_variable$rho.c[, 2], digits = 2), " - ", 
               round(ccc_variable$rho.c[, 3], digits = 2), "]", sep = "")
  
  # Create the PNG file name based on the combination
  png_filename <- paste0("concordance_correlation_plot_Inter_", comb, ".png")
  png(png_filename, width = 200, height = 200, units = "mm", res = 1000)
  
  # Plotting
  plot(var1, var2, 
       xlab = paste("Set", sets[1]), 
       ylab = paste("Set", sets[2]), 
       pch = 16,
       main = paste("Inter-observer comparison: Set", sets[1], "and Set", sets[2]),
       cex.main = 1.6, cex.lab = 1.3, cex.axis = 0.8)
  
  # Add line of perfect concordance
  abline(a = 0, b = 1, lty = 2, col = "red", lwd = 2)
  
  # Calculate RMA slope and intercept
  sd_x <- sd(var1)
  sd_y <- sd(var2)
  rma_slope <- sd_y / sd_x
  rma_intercept <- mean(var2) - rma_slope * mean(var1)
  
  # Add the RMA line
  abline(a = rma_intercept, b = rma_slope, lty = 1, lwd = 2)
  
  # Add legend
  legend("topleft", legend = c("Line of perfect concordance", "Reduced major axis"), 
         lty = c(2, 1), col = c("red", "black"), lwd = c(2, 2), bty = "n",
         cex = 1)
  
  # Position text in the lower right corner
  x_pos <- max(var1) * 0.8
  y_pos <- min(var2)
  
  # Add text
  text(x = x_pos, y = y_pos, labels = lab, pos = 4, cex = 0.9, col = "black")

  dev.off()
}


# 3D Area: Export statistical results -------------------------------------

### Automatically export the heatmap as a PDF file ###
ggsave("pairwise_ccc_heatmap.pdf", plot = heatmap_plot, width = 200, height = 170, units = "mm", dpi = 1000)


### Export intra-observer statistics as excel file ###
# Create a new workbook
wb <- createWorkbook()

# Write data frames to worksheets
addWorksheet(wb, "Descriptive Statistics")
writeData(wb, "Descriptive Statistics", area_intra_descr)

addWorksheet(wb, "Percentage error")
writeData(wb, "Percentage error", area_intra_perc)

# Save the workbook
saveWorkbook(wb, file = file.path(getwd(), "3DAreas_intra.xlsx"), overwrite = TRUE)


### Export inter_observer results as excel file ###
# Create a new workbook
wb <- createWorkbook()

# Write data frames to worksheets
addWorksheet(wb, "Descriptive statistics")
writeData(wb, "Descriptive statistics", area_inter_descr)

addWorksheet(wb, "Percentage error")
writeData(wb, "Percentage error", area_inter_perc)

addWorksheet(wb, "Percentage error with average")
writeData(wb, "Percentage error with average", area_inter_perc_av)

# Save the workbook
saveWorkbook(wb, file = file.path(getwd(), "3DAreas_inter.xlsx"), overwrite = TRUE)


### Export correlation coefficient results as excel file ###
# Create a new workbook
wb <- createWorkbook()

# Write data frames to worksheets
addWorksheet(wb, "Area normality")
writeData(wb, "Area normality", area_normality)

addWorksheet(wb, "ICC")
writeData(wb, "ICC", area_icc_result)

addWorksheet(wb, "Lin's CCC")
writeData(wb, "Lin's CCC", area_ccc_result)

addWorksheet(wb, "Lin's CCC correlation matrix")
writeData(wb, "Lin's CCC correlation matrix", ccc_matrix)

# Save the workbook
saveWorkbook(wb, file = file.path(getwd(), "3DAreas_CorrelationCoefficient.xlsx"), overwrite = TRUE)

