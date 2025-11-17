# Calculating the Alpha and Beta Diversity for the Metagenomic Samples
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1


#set your working directory were everything will be stored and saved
setwd("C:/Users/Sandersonh/OneDrive - AGR-AGR/Desktop/Wen Chen/Project_1_landuse_AMR_virulence/Statistics_Analysis")


#install.packages("tidyverse")
#install.packages("vegan")
#install.packages(c("dplyr"))

library(dplyr)
library(tidyverse)
library(vegan)

#getting genus level relative abundances
genus_data_nn <- read.delim("genus.tsv", header = TRUE, stringsAsFactors = FALSE)
genus_data_relabund_nn <- genus_data_nn %>%
  group_by(SAMPLE_ID) %>%
  mutate(RelativeAbundance = READS_COUNT / sum(READS_COUNT))

cross_table_genus_nn <- genus_data_relabund_nn %>%
  select(SAMPLE_ID, TAXONOMY, RelativeAbundance) %>%
  spread(TAXONOMY, RelativeAbundance, fill = 0)


#Removing SAmples we didn't use
sample_ids_to_remove <- c(
  "SN253-20210517", "SN253-20210531", "SN253-20210628", "SN253-20210712",
  "SN253-20211004", "SN253-20211018", "SN253-20170704", "SN253-20170717",
  "SN253-20171023", "SN253-20171106", "SN253-20180430", "SN253-20180528",
  "SN253-20180807", "SN253-20180820", "SN253-20181113", "SN253-20190429",
  "SN253-20190513", "SN253-20190527", "SN253-20190624", "SN253-20190708",
  "SN253-20190722", "SN253-20190812", "SN253-20190923", "SN253-20191021",
  "SN253-20191104", "SN253-20201005", "SN253-20201019", "SN253-20201102",
  "SN253-20170523", "SN253-20181015", "SN005-20201007", "SN006-20201007",
  "SN010-20201007", "SN018-20180418", "SN018-20201007", "SN019-20180418",
  "SN019-20201007", "SN020-20180418", "SN020-20201007", "SN021-20180418",
  "SN021-20201007", "SN024-20201007", "SN20-20201019", "SN20-20201102",
  "SN253-20170605", "SN253-20201007"
)

filtered_table_genus <- cross_table_genus_nn %>%
  filter(!SAMPLE_ID %in% sample_ids_to_remove)


#filter sample ID based on another dataframe with the right number of sample_IDs

data_sampleID <- read.csv("table_lme_metadata_AMRVF_taxa_month.csv", header = TRUE, stringsAsFactors = FALSE)

filtered_table_genus <- filtered_table_genus %>%
  filter(SAMPLE_ID %in% data_sampleID$SAMPLE_ID)




colnames(filtered_table_genus)[colnames(filtered_table_genus) == "Pusillimonas_(ex Stolz et al. 2005)"] <- "Pusillimonas"

#Get sample site IDs from sample_ids
filtered_table_genus$Sample_Site <- sub("-(\\d+)$", "", filtered_table_genus$SAMPLE_ID)


#Create the Blocks
filtered_table_genus$Block <- ifelse(filtered_table_genus$Sample_Site %in% c("SN018", "SN019"), "Block1819",
                               ifelse(filtered_table_genus$Sample_Site %in% c("SN020", "SN021"), "Block2021", paste0("Block_", filtered_table_genus$Sample_Site)))

#remove SAmple_ID and Sample_Site columns and make Block the row.names
# Remove Sample_ID and Sample_Site columns
filtered_table_genus <- filtered_table_genus[, !(names(filtered_table_genus) %in% c("Sample_Site"))]

# Set Block column as row names
rownames(filtered_table_genus) <- filtered_table_genus$SAMPLE_ID
filtered_table_genus$SAMPLE_ID <- NULL


#Alpha Diversity
#remove Block column to do alphadiversity

# Select all columns except the last one
selected_columns <- filtered_table_genus[, -ncol(filtered_table_genus)]

# Assuming 'Sample_ID' is the last column in your dataframe, you can use the selected columns in the diversity function
alpha_diversity_results_shannon <- diversity(selected_columns, index = "shannon")
alpha_diversity_results_simpson <- diversity(selected_columns, index = "simpson")

# Extract evenness and richness values
evenness <- alpha_diversity_results_shannon 
richness <- alpha_diversity_results_simpson
#Will need to figure this out for by Block and Land_use and compare alpha diversity
# Create a data frame for plotting
metadata <- read.csv("Metadata_clean_AMRVFMG_project.csv", header = TRUE, stringsAsFactors = FALSE) 

metadata <- metadata %>%
  rename(SAMPLE_ID = Sample_ID_Genome)

filtered_metadata <- metadata %>%
  filter(SAMPLE_ID %in% data_sampleID$SAMPLE_ID)


diversity_df <- data.frame(
  SAMPLE_ID = filtered_metadata$SAMPLE_ID,
  Block = filtered_table_genus$Block, 
  Land_Use = filtered_metadata$Land_use,
  Evenness = evenness,   
  Richness = richness     
)




# Kruskal-Wallis test
kruskal_results_evenness <- kruskal.test(Evenness ~ Block, data = diversity_df)
kruskal_results_richness <- kruskal.test(Richness ~ Block, data = diversity_df)

# Post-hoc tests for pairwise comparisons (Wilcoxon signed-rank test)
posthoc_evenness <- pairwise.wilcox.test(diversity_df$Evenness, diversity_df$Sample, p.adj = "bonferroni")
posthoc_richness <- pairwise.wilcox.test(diversity_df$Richness, diversity_df$Sample, p.adj = "bonferroni")

# Create bar plots using ggplot2 with annotations for significant differences
#install.packages("ggpubr")
#install.packages("ggsignif")
library(ggplot2)
library(ggpubr)
library(ggsignif)


# Evenness Block
# Create a box and whisker plot for Evenness
# Create a box and whisker plot for Evenness

# Create a box and whisker plot for Evenness
evenness_boxplot <- ggplot(diversity_df, aes(x = Block, y = Evenness)) +
  geom_boxplot() +
  labs(title = "Evenness-Block Comparison", x = "Block", y = "Evenness")
my_comparisons <- list( c("Block_SN005", "Block1819"), c("Block_SN006", "Block1819"), c("Block_SN006", "Block2021"), c("Block_SN006", "Block_SN024"), c("Block_SN010", "Block1819"))
evenness_boxplot + stat_compare_means(comparisons = my_comparisons)

#Richness Block
richness_boxplot <- ggplot(diversity_df, aes(x = Block, y = Richness)) +
  geom_boxplot() +
  labs(title = "Evenness-Block Comparison", x = "Block", y = "Richness")
my_comparisons <- list(c("Block_SN005", "Block_SN006"), c("Block_SN006", "Block1819"))
richness_boxplot + stat_compare_means(comparisons = my_comparisons)  



#Evenness Land use
# Create a box and whisker plot for Evenness
diversity_df$Land_Use <- as.factor(diversity_df$Land_Use)


evenness_boxplot <- ggplot(diversity_df, aes(x = Land_Use, y = Evenness)) +
  geom_boxplot() +
  labs(title = "Evenness-Land_Use Comparison", x = "Land_Use", y = "Evenness")
my_comparisons <- list(c("Agri_ditch", "Mixed"), c("Agri_ditch", "Forested"), c("Mixed", "Forested"))
# Add pairwise comparisons
evenness_boxplot + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", grouping.var = "Land_Use")

#Richness Land use
richness_boxplot <- ggplot(diversity_df, aes(x = Land_Use, y = Richness)) +
  geom_boxplot() +
  labs(title = "Richness-Land_Use Comparison", x = "Land_Use", y = "Richness")
my_comparisons <- list(c("Agri_ditch", "Mixed"), c("Agri_ditch", "Forested"), c("Mixed", "Forested"))
# Add pairwise comparisons
richness_boxplot + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", grouping.var = "Land_Use")



# Print the plots
print(evenness_boxplot)
print(richness_boxplot)

#Getting pvalues in atable with the pairs
#Evenness_Block
# Get unique blocks
unique_blocks <- unique(diversity_df$Block)

# Initialize an empty dataframe for all pairs
all_pairs <- data.frame(Block1 = character(), Block2 = character(), p.value = numeric(), stringsAsFactors = FALSE)

# Perform pairwise Wilcoxon rank sum test for each pair of blocks
for (block1 in unique_blocks) {
  for (block2 in unique_blocks) {
    if (block1 < block2) {
      values_block1 <- diversity_df$Evenness[diversity_df$Block == block1]
      values_block2 <- diversity_df$Evenness[diversity_df$Block == block2]
      
      # Check if both blocks have enough observations
      if (length(values_block1) > 1 && length(values_block2) > 1) {
        wilcox_result <- wilcox.test(values_block1, values_block2)
        p_val <- wilcox_result$p.value
        all_pairs <- rbind(all_pairs, data.frame(Block1 = block1, Block2 = block2, p.value = p_val))
      }
    }
  }
}

# Print or use 'all_pairs' dataframe as needed
print(all_pairs)
Evenness_Block_pvalues <- all_pairs

#Richness_Block
# Get unique blocks
unique_blocks <- unique(diversity_df$Block)

# Initialize an empty dataframe for all pairs
all_pairs <- data.frame(Block1 = character(), Block2 = character(), p.value = numeric(), stringsAsFactors = FALSE)

# Perform pairwise Wilcoxon rank sum test for each pair of blocks
for (block1 in unique_blocks) {
  for (block2 in unique_blocks) {
    if (block1 < block2) {
      values_block1 <- diversity_df$Richness[diversity_df$Block == block1]
      values_block2 <- diversity_df$Richness[diversity_df$Block == block2]
      
      # Check if both blocks have enough observations
      if (length(values_block1) > 1 && length(values_block2) > 1) {
        wilcox_result <- wilcox.test(values_block1, values_block2)
        p_val <- wilcox_result$p.value
        all_pairs <- rbind(all_pairs, data.frame(Block1 = block1, Block2 = block2, p.value = p_val))
      }
    }
  }
}

# Print or use 'all_pairs' dataframe as needed
print(all_pairs)
Richness_Block_pvalues <- all_pairs

#Evenness Landuse

# Get unique blocks
unique_blocks <- unique(diversity_df$Land_Use)

# Initialize an empty dataframe for all pairs
all_pairs <- data.frame(Block1 = character(), Block2 = character(), p.value = numeric(), stringsAsFactors = FALSE)

# Perform pairwise Wilcoxon rank sum test for each pair of blocks
for (block1 in unique_blocks) {
  for (block2 in unique_blocks) {
    if (block1 < block2) {
      values_block1 <- diversity_df$Evenness[diversity_df$Land_Use == block1]
      values_block2 <- diversity_df$Evenness[diversity_df$Land_Use == block2]
      
      # Check if both blocks have enough observations
      if (length(values_block1) > 1 && length(values_block2) > 1) {
        wilcox_result <- wilcox.test(values_block1, values_block2)
        p_val <- wilcox_result$p.value
        all_pairs <- rbind(all_pairs, data.frame(Block1 = block1, Block2 = block2, p.value = p_val))
      }
    }
  }
}

# Print or use 'all_pairs' dataframe as needed
print(all_pairs)
Evenness_Landuse__pvalues <- all_pairs

#Evenness Landuse

# Get unique blocks
unique_blocks <- unique(diversity_df$Land_Use)

# Initialize an empty dataframe for all pairs
all_pairs <- data.frame(Block1 = character(), Block2 = character(), p.value = numeric(), stringsAsFactors = FALSE)

# Perform pairwise Wilcoxon rank sum test for each pair of blocks
for (block1 in unique_blocks) {
  for (block2 in unique_blocks) {
    if (block1 < block2) {
      values_block1 <- diversity_df$Richness[diversity_df$Land_Use == block1]
      values_block2 <- diversity_df$Richness[diversity_df$Land_Use == block2]
      
      # Check if both blocks have enough observations
      if (length(values_block1) > 1 && length(values_block2) > 1) {
        wilcox_result <- wilcox.test(values_block1, values_block2)
        p_val <- wilcox_result$p.value
        all_pairs <- rbind(all_pairs, data.frame(Block1 = block1, Block2 = block2, p.value = p_val))
      }
    }
  }
}

# Print or use 'all_pairs' dataframe as needed
print(all_pairs)
Richness_Landuse__pvalues <- all_pairs





#Beta Diversity
#install.packages("compositions")
library(compositions)

taxa_data_genus <- selected_columns

install.packages(c("vegan", "GUniFrac", "ggplot2", "ellipse"))
library(vegan)
library(GUniFrac)
library(ggplot2)
library(ellipse)
clr_transformed_data <- clr(taxa_data_genus)
# Calculate Bray-Curtis dissimilarity matrix
euclidean_dist <- vegdist(clr_transformed_data, method = "euclidean")

# Perform PCoA
pcoa_result <- cmdscale(bray_curtis_dist)


# Create an ordination plot with ellipses
ordination_df <- data.frame(
  SampleID = sample_ids,
  PC1 = pcoa_result[, 1],
  PC2 = pcoa_result[, 2],
  LandUse = filtered_metadata$Land_use,
  Block = filtered_table_genus$Block
)

ggplot(ordination_df, aes(x = PC1, y = PC2, color = LandUse, shape = Block)) +
  geom_point() +
  stat_ellipse(geom = "polygon", aes(fill = LandUse), alpha = 0.05) +
  labs(title = "PCoA - Euclidean Dissimilarity", color = "Land Use", shape = "Block") +
  theme_minimal()

# PERMANOVA between land use classes
# Perform PERMANOVA
# Perform PERMANOVA
pc_data <- ordination_df[, c("PC1", "PC2")]
permanova_result <- adonis2(pc_data ~ LandUse + Block, data = ordination_df, method = "euclidean")

# Print the results
print(permanova_result)

land_use_classes <- c("Agri_ditch", "Forested", "Mixed")


# Subset the data for the specified land use classes
# Create a distance matrix
dist_matrix <- dist(subset_data[, c("PC1", "PC2")])

# Perform pairwise comparisons between the specified land use classes manually
pairwise_results <- list()
for (i in 1:(length(land_use_classes) - 1)) {
  for (j in (i + 1):length(land_use_classes)) {
    group1 <- land_use_classes[i]
    group2 <- land_use_classes[j]
    
    # Subset data for the current pair
    subset_pair <- subset_data[subset_data$LandUse %in% c(group1, group2), ]
    
    # Create a new distance matrix for the current pair
    dist_matrix_pair <- dist(subset_pair[, c("PC1", "PC2")])
    
    # Perform PERMANOVA using adonis
    pairwise_results[[paste(group1, group2, sep = "_vs_")]] <- adonis2(dist_matrix_pair ~ LandUse, data = subset_pair, method = "euclidean")
  }
}

# Print the results
for (i in seq_along(pairwise_results)) {
  cat("Pairwise Comparison:", names(pairwise_results)[i], "\n")
  print(pairwise_results[[i]])
  cat("\n")
}

#permanova between blocks
block_classes <- c("Block_SN005", "Block_SN006", "Block_SN010", "Block_SN024", "Block1819", "Block2021")


# Subset the data for the specified land use classes
# Create a distance matrix
dist_matrix <- dist(subset_data[, c("PC1", "PC2")])

# Perform pairwise comparisons between the specified land use classes manually
pairwise_results <- list()
for (i in 1:(length(block_classes) - 1)) {
  for (j in (i + 1):length(block_classes)) {
    group1 <- block_classes[i]
    group2 <- block_classes[j]
    
    # Subset data for the current pair
    subset_pair <- subset_data[subset_data$Block %in% c(group1, group2), ]
    
    # Create a new distance matrix for the current pair
    dist_matrix_pair <- dist(subset_pair[, c("PC1", "PC2")])
    
    # Perform PERMANOVA using adonis
    pairwise_results[[paste(group1, group2, sep = "_vs_")]] <- adonis2(dist_matrix_pair ~ Block, data = subset_pair, method = "euclidean")
  }
}

# Print the results
for (i in seq_along(pairwise_results)) {
  cat("Pairwise Comparison:", names(pairwise_results)[i], "\n")
  print(pairwise_results[[i]])
  cat("\n")
}

