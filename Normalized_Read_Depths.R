# GEtting Normalized REad Depths for each metagenomic sample and analyzing it
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1
#set your working directory were everything will be stored and saved
setwd("C:/Users/Sandersonh/OneDrive - AGR-AGR/Desktop/Wen Chen/Project_1_landuse_AMR_virulence/Statistics_Analysis")

readcounts <- read.csv("Metagenomic_Samples_Read_Counts.csv", header = TRUE, stringsAsFactors = FALSE)

#get normalized read depths from Read Counts
# Define function for RPM normalization
normalize_rpm <- function(read_counts) {
  total_reads <- sum(read_counts)
  normalized_depths <- (read_counts / total_reads) * 1e6
  return(normalized_depths)
}

readcounts$Normalized_Read_Depths <- normalize_rpm(readcounts$Total_Read_Count)

#Read depths per sample block

metadata <- read.csv("Metadata_clean_AMRVFMG_project.csv", header = TRUE, stringsAsFactors = FALSE)


library(tidyr)
library(dplyr)

selected_columns <- c("Sample_ID_Genome", "Sample_Site")
metadata_dataframe <- metadata[selected_columns]
metadata_dataframe <- metadata_dataframe %>% rename(SAMPLE_ID = Sample_ID_Genome)
combined_df <- merge(metadata_dataframe, readcounts, by = "SAMPLE_ID", all = TRUE)
combined_df <- na.omit(combined_df)

combined_df$Block <- ifelse(combined_df$Sample_Site %in% c("SN018", "SN019"), "Block1819",
                            ifelse(combined_df$Sample_Site %in% c("SN020", "SN021"), "Block2021", paste0("Block_", combined_df$Sample_Site)))
library(ggplot2)
library(ggpubr)
#make a box plot out of the read depths for each block
p= ggplot(data = combined_df, aes(x = Block, y = Normalized_Read_Depths)) +
  geom_boxplot() +
  labs(x = "Block", y = "Normalized Read Depth") +
  ggtitle("Box Plot of Normalized Read Depth by Block")


#see if there is a significant difference between read depth between blocks
pairwise_pvalues <- list()
blocks <- unique(combined_df$Block)

for (i in 1:(length(blocks)-1)) {
  for (j in (i+1):length(blocks)) {
    block1 <- blocks[i]
    block2 <- blocks[j]
    subset_data <- combined_df[combined_df$Block %in% c(block1, block2),]
    p_val <- wilcox.test(Normalized_Read_Depths ~ Block, data = subset_data)$p.value
    comparison <- paste(block1, block2, sep = "_vs_")
    pairwise_pvalues[[comparison]] <- p_val
  }
}

# Get min, max and average normalized read Depth
highest_readdepth <- max(combined_df$Normalized_Read_Depths)
lowest_readdepth <- min(combined_df$Normalized_Read_Depths)
average_readdepth <- mean(combined_df$Normalized_Read_Depths)

#Print all the results

print(p)
print(pairwise_pvalues)
print(highest_readdepth)
print(lowest_readdepth)
print(average_readdepth)
