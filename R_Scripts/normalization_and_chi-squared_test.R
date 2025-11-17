# Draft script for normalizing and transforming metagenomic taxonomic data
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1

#set your working directory were everything will be stored and saved
setwd("C:/Users/Sandersonh/OneDrive - AGR-AGR/Desktop/Wen Chen/Project_1_landuse_AMR_virulence/Statistics_Analysis")


#if(!requireNamespace("BiocManager")){
#  install.packages("BiocManager")
#}
#BiocManager::install("phyloseq")

# Read data from your TSV file
genus_data <- read.delim("genus_normalized.tsv", header = TRUE, stringsAsFactors = FALSE)

genus_data_relabund <- genus_data %>%
  group_by(SAMPLE_ID) %>%
  mutate(RelativeAbundance = READS_COUNT / sum(READS_COUNT) * 100)


species_data <- read.delim("species_normalized.tsv", header = TRUE, stringsAsFactors = FALSE)

species_data_relabund <- species_data %>%
  group_by(SAMPLE_ID) %>%
  mutate(RelativeAbundance = READS_COUNT / sum(READS_COUNT) * 100)

family_data <- read.delim("family_normalized.tsv", header = TRUE, stringsAsFactors = FALSE)
family_data_relabund <- family_data %>%
  group_by(SAMPLE_ID) %>%
  mutate(RelativeAbundance = READS_COUNT / sum(READS_COUNT) * 100)



order_data <- read.delim("order_normalized.tsv", header = TRUE, stringsAsFactors = FALSE)
order_data_relabund <- order_data %>%
  group_by(SAMPLE_ID) %>%
  mutate(RelativeAbundance = READS_COUNT / sum(READS_COUNT) * 100)



#install.packages("tidyr")
#install.packages("dplyr")


library(tidyr)
library(dplyr)

cross_table_species <- species_data_relabund %>%
  select(SAMPLE_ID, TAXONOMY, RelativeAbundance) %>%
  spread(TAXONOMY, RelativeAbundance, fill = 0)

cross_table_genus <- genus_data_relabund %>%
  select(SAMPLE_ID, TAXONOMY, RelativeAbundance) %>%
  spread(TAXONOMY, RelativeAbundance, fill = 0)

cross_table_family <- family_data_relabund %>%
  select(SAMPLE_ID, TAXONOMY, RelativeAbundance) %>%
  spread(TAXONOMY, RelativeAbundance, fill = 0)

cross_table_order <- order_data_relabund %>%
  select(SAMPLE_ID, TAXONOMY, RelativeAbundance) %>%
  spread(TAXONOMY, RelativeAbundance, fill = 0)


write.csv(cross_table_species, file = "species_cross_table.csv", row.names = FALSE)
write.csv(cross_table_genus, file = "genus_cross_table.csv", row.names = FALSE)
write.csv(cross_table_family, file = "family_cross_table.csv", row.names = FALSE)
write.csv(cross_table_order, file = "order_cross_table.csv", row.names = FALSE)

# create one big table with all the taxonomic data

merged_table <- merge(cross_table_order, cross_table_family, by = "SAMPLE_ID", all = TRUE)
merged_table <- merge(merged_table, cross_table_genus, by = "SAMPLE_ID", all = TRUE)


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

filtered_table <- merged_table %>%
  filter(!SAMPLE_ID %in% sample_ids_to_remove)


# Save the merged table as a CSV file
write.csv(merged_table, file = "taxonomy_merged_table.csv", row.names = FALSE)
# Save the filtered table as a CSV file
write.csv(filtered_table, file = "taxonomy_filtered_table.csv", row.names = FALSE)

#Import AMR and VF data
AMR <- read.csv("Read_mapping_AMR_crosstab.csv", header = TRUE, stringsAsFactors = FALSE) 
VF <- read.csv("Read_mapping_VFDB_crosstab.csv", header = TRUE, stringsAsFactors = FALSE)

#merge AMR and VF tables
merged_table_AMRVF <- merge(AMR, VF, by = "Samples", all = TRUE)
#Change samples column to SAMPLE_ID
merged_table_AMRVF <- merged_table_AMRVF %>%
  rename(SAMPLE_ID = Samples)
#filter out sample_IDs that were removed from metadata
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

filtered_table_AMR_VF <- merged_table_AMRVF %>%
  filter(!SAMPLE_ID %in% sample_ids_to_remove)

write.csv(filtered_table_AMR_VF, file = "filtered_AMRVF_crosstab.csv", row.names = FALSE)

#Import metadata table
metadata <- read.csv("Metadata_clean_AMRVFMG_project.csv", header = TRUE, stringsAsFactors = FALSE) 

metadata <- metadata %>%
  rename(SAMPLE_ID = Sample_ID_Genome)


#merge taxonomy, AMRVF and metadata tables together
merged_table_all <- merge(metadata, filtered_table_AMR_VF, by = "SAMPLE_ID", all = TRUE)
merged_table_all <- merge(merged_table_all, filtered_table, by = "SAMPLE_ID", all = TRUE)

#save complete overall table
write.csv(filtered_table_AMR_VF, file = "merged_overall_table_metadata_AMRVF_taxonomy.csv", row.names = FALSE)

#install.packages("rstatix")
#install.packages("purrr")
library("rstatix")
library("purrr")


##Pearson Correlations: land use vs AMR vs VF vs taxa
filtered_table_AMRVF_column_names <- colnames(filtered_table_AMR_VF)
print(filtered_table_AMRVF_column_names)
# Display all column names for taxonomy_filtered_table
taxonomy_filtered_table_column_names <- colnames(filtered_table)
print(taxonomy_filtered_table_column_names)
# Select all columns except SAMPLE_ID and other non-desirable columns from the tables
selected_columns_AMRVF <- setdiff(filtered_table_AMRVF_column_names, "SAMPLE_ID")
selected_columns_taxonomy <- setdiff(taxonomy_filtered_table_column_names, "SAMPLE_ID")


selected_columns_AMRVF <-setdiff(selected_columns_AMRVF, "Total_AMR_per_Sample") 
selected_columns_AMRVF <-setdiff(selected_columns_AMRVF, "Total_per_Samples")

# Combine selected columns from both tables
selected_columns <- c("Land_use", selected_columns_taxonomy)

# Select relevant columns for correlation analysis
cor_data <- merged_table_all %>%
  select(all_of(selected_columns))

# Filter only numeric columns for Kruskal-Wallis test
numeric_columns <- cor_data %>%
  select(where(is.numeric))

# Combine Land_use with numeric columns

data_for_kruskal <- numeric_columns %>%
  gather(key = "Variable", value = "Value") %>%
  cbind(Land_use = cor_data$Land_use)

# Convert Land_use to a factor
data_for_kruskal$Land_use <- as.factor(data_for_kruskal$Land_use)

# Kruskal-Wallis test
kruskal_results <- kruskal_test(Value ~ Land_use, data = data_for_kruskal)

# Print the Kruskal-Wallis test results
print(kruskal_results)
# Post-hoc pairwise comparisons
posthoc_results <- data_for_kruskal %>%
  pairwise_wilcox_test(
    Value ~ Land_use,
    p.adjust.method = "BH"
  )
#save posthoc results as Dataframe
posthoc_results_df <- as.data.frame(posthoc_results)

# Save as a CSV file
write.csv(posthoc_results_df, file = "TAXA_posthoc_results.csv", row.names = FALSE)

#land-use vs AMR and VF
selected_data <- merged_table_all %>%
  select("Land_use", all_of(selected_columns_AMRVF))

chi_squared_results <- list()

# Perform chi-squared test for each selected column
for (column_name in selected_columns_AMRVF) {
  contingency_table <- table(selected_data$Land_use, selected_data[[column_name]])
  chi_squared_test <- chisq.test(contingency_table)
  chi_squared_results[[column_name]] <- chi_squared_test
}

# Display the chi-squared test results
print(chi_squared_results)

# Create an empty data frame to store results
chisq_results_df <- data.frame(Column = character(), 
                               ChiSquaredStatistic = numeric(),
                               PValue = numeric(),
                               stringsAsFactors = FALSE)

# Extract information from each chi-squared test
for (column_name in selected_columns_AMRVF) {
  chi_squared_test <- chi_squared_results[[column_name]]
  
  # Extract relevant information
  result_row <- data.frame(
    Column = column_name,
    ChiSquaredStatistic = chi_squared_test$statistic,
    PValue = chi_squared_test$p.value
  )
  
  # Append to the results data frame
  chisq_results_df <- bind_rows(chisq_results_df, result_row)
}

# Save the data frame as a CSV file
write.csv(chisq_results_df, file = "AMR_VF_Landuse_chi_squared_results.csv", row.names = FALSE)

# Display the data frame
print(chisq_results_df)


#####Season vs AMR/VF#################
#land-use vs AMR and VF
selected_data <- merged_table_all %>%
  select("Season", all_of(selected_columns_AMRVF))

chi_squared_results <- list()

# Perform chi-squared test for each selected column
for (column_name in selected_columns_AMRVF) {
  contingency_table <- table(selected_data$Season, selected_data[[column_name]])
  chi_squared_test <- chisq.test(contingency_table)
  chi_squared_results[[column_name]] <- chi_squared_test
}

# Display the chi-squared test results
print(chi_squared_results)

# Create an empty data frame to store results
chisq_results_df <- data.frame(Column = character(), 
                               ChiSquaredStatistic = numeric(),
                               PValue = numeric(),
                               stringsAsFactors = FALSE)

# Extract information from each chi-squared test
for (column_name in selected_columns_AMRVF) {
  chi_squared_test <- chi_squared_results[[column_name]]
  
  # Extract relevant information
  result_row <- data.frame(
    Column = column_name,
    ChiSquaredStatistic = chi_squared_test$statistic,
    PValue = chi_squared_test$p.value
  )
  
  # Append to the results data frame
  chisq_results_df <- bind_rows(chisq_results_df, result_row)
}

# Save the data frame as a CSV file
write.csv(chisq_results_df, file = "AMR_VF_Season_chi_squared_results.csv", row.names = FALSE)

# Display the data frame
print(chisq_results_df)

#Differential Analysis of Taxa between Samples




#Visualizing taxonomic data (not using)

# Load the ggplot2 package
#if (!requireNamespace("ggplot2", quietly = TRUE)) {
 # install.packages("ggplot2")
#}
library(ggplot2)

# Create a stacked barplot
ggplot(genus_data_relabund, aes(x = SAMPLE_ID, y = RelativeAbundance, fill = TAXONOMY)) +
  geom_bar(stat = "identity") +
  labs(title = "Stacked Barplot of Genus Relative Abundance",
       x = "Sample ID", y = "Relative Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#filter RA by a cutoff
filtered_RA_data <- genus_data_relabund %>%
  filter(RelativeAbundance >= 0.1)

# Create a stacked barplot with the filtered data
plot <- ggplot(filtered_RA_data, aes(x = SAMPLE_ID, y = RelativeAbundance, fill = TAXONOMY)) +
  geom_bar(stat = "identity") +
  labs(title = "Stacked Barplot of Genus Relative Abundance (>= 0.1%)",
       x = "Sample ID", y = "Relative Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",  # Adjust the legend position if needed
        legend.text = element_text(size = 10 / 10))  # Reduce the size of the legend by a factor of 10

# Display the plot
ggsave("Genus_stacked_barplot.jpg", plot, width = 10, height = 6, units = "in")
