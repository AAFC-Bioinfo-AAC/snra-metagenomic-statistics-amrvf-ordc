# Eigengene Network Analysis
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1

#set your working directory were everything will be stored and saved
setwd("C:/Users/Sandersonh/OneDrive - AGR-AGR/Desktop/Wen Chen/Project_1_landuse_AMR_virulence/Statistics_Analysis")


#install.packages('WGCNA')
#install.packages('matrixStats')
#install.packages('impute')
#install.packages("RSQLite")
#install.packages("BiocManager")
#BiocManager::install("GO.db")
#BiocManager::install("WGCNA")
library(RSQLite)
library(matrixStats)
library(impute)
library(BiocManager)
library(WGCNA)

# Load and preprocess data
data_file <- "table_lme_metadata_AMRVF_taxa_month.csv"

# Load your metagenomic taxonomic data
taxonomic_data <- read.csv(data_file, header = TRUE, row.names = 1)

AMR_VF_column_names <- colnames(AMR_VF)
order_column_names <- colnames(clr_data_order_df)
family_column_names <- colnames(clr_data_family_df)
genus_column_names <- colnames(clr_data_genus_df)
selected_columns <- c("Land_use", "Season")
metadata_dataframe <- metadata[selected_columns]
metadata_column_names <- colnames(metadata_dataframe)  

exclude_columns <- c(order_column_names, family_column_names)
include_columns <- c(metadata_column_names, AMR_VF_column_names, genus_column_names)

# Filter taxonomic data based on include and exclude columns
filtered_taxonomic_data <- taxonomic_data[, colnames(taxonomic_data) %in% include_columns & !colnames(taxonomic_data) %in% exclude_columns]

# Display the filtered taxonomic data
print(filtered_taxonomic_data)


#Transform and Filter data
# Make sure Land_use, Season and AMRVF data is seen as factors
factor_columns = c("Land_use", "Season", AMR_VF_column_names)
factor_columns <- colnames(filtered_taxonomic_data) %in% factor_columns
filtered_taxonomic_data[, factor_columns] <- lapply(filtered_taxonomic_data[, factor_columns], as.factor)
#Filter unwanted Columns for AMRVF dataset that are still in the filtered dataset
excludeAMRcolumns = c("Total_AMR_per_Sample", "Total_per_Samples")
filtered_taxonomic_data <- filtered_taxonomic_data[, !(colnames(filtered_taxonomic_data) %in% excludeAMRcolumns)]

#Log-transform the non-factor data is not required because the data is already clr transformed

#transpose dataset
filtered_taxonomic_data <- t(filtered_taxonomic_data)

#Make the dataframe that is only numeric
exclude_columns<- c("Land_use", "Season", AMR_VF_column_names)
filtered_taxonomic_data_numeric <- filtered_taxonomic_data[,!(colnames(filtered_taxonomic_data) %in% exclude_columns)]
filtered_taxonomic_data_numeric <- as.matrix(filtered_taxonomic_data_numeric)
# Convert to numeric if possible, replacing non-numeric values with NA
filtered_taxonomic_data_numeric <- apply(filtered_taxonomic_data_numeric, 2, function(x) as.numeric(as.character(x)))

# Check if 'x' is numeric now
is.numeric(filtered_taxonomic_data_numeric)


library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(magrittr)      # provides the %>% operator
library(WGCNA)
#For sample network
#getting the right orientation of the dataset for samples to make the network
filtered_taxonomic_data_sample <- t(filtered_taxonomic_data)
exclude_columns<- c("Land_use", "Season", AMR_VF_column_names)
# Find common elements between colnames and exclude_columns
common_columns <- intersect(colnames(filtered_taxonomic_data_sample), exclude_columns)
# Exclude columns with common names
filtered_taxonomic_data_numeric_sample <- filtered_taxonomic_data_sample[ , !colnames(filtered_taxonomic_data_sample) %in% common_columns]
filtered_taxonomic_data_numeric_sample <- as.matrix(filtered_taxonomic_data_numeric_sample)
# Convert to numeric if possible, replacing non-numeric values with NA
filtered_taxonomic_data_numeric_sample <- apply(filtered_taxonomic_data_numeric_sample, 2, function(x) as.numeric(as.character(x)))
#Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
#call the network topology analysis function
sft= pickSoftThreshold(filtered_taxonomic_data_numeric_sample, powerVector = powers, verbose = 5)

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")


picked_power = 8
temp_cor <- cor
cor <- WGCNA::cor

#iltered_taxonomic_data_numeric <- t(filtered_taxonomic_data_numeric)
netwk_samples<- blockwiseModules(filtered_taxonomic_data_numeric_sample, power = picked_power, networkType = "signed", deepSplit = 4, pamRespectsDendro = F, minModuleSize =2, maxBlockSize = 4000, reassignThreshold = 0, mergeCutHeight = 0.05, saveTOMs = T, saveTOMFileBase = "ER", numericLabels = T, verbose = 3)

mergedColors = labels2colors(netwk_samples$colors)

png("EN_dendro_colours.png", width = 800, height = 600)
plotDendroAndColors(netwk_samples$dendrograms[[1]], mergedColors[netwk_samples$blockGenes[[1]]], "Module Colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()


#Get module Eigengenes per cluster (Genus)

module_df <- data.frame(samples = names(netwk_samples$colors), colors = labels2colors(netwk_samples$colors))
module_df[1:5,]
#seems to be genus taxa?
write_delim(module_df, file = "genus_modules.txt", delim = "\t")


#GEt Module Eigengenes per cluster (samples)

MEs0 <- moduleEigengenes(filtered_taxonomic_data_numeric_sample, mergedColors)$eigengenes
#reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)
#add sample names
MEs0$sample = row.names(t(filtered_taxonomic_data))
#tidy & plot data
mME = MEs0 %>%
  pivot_longer(-sample) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )
png("EN_Module_Sample_Relationship_heatmap.png", width = 800, height = 600)
mME %>% ggplot(., aes(x=sample, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, size=6)) +
  labs(title = "Module-Sample Relationships", y = "Modules", fill="corr")
dev.off()

#making a dendrogram with land_use and season for samples network
df_labels <- taxonomic_data[, c("Land_use", "Season")]
install.packages("pheatmap")

library(pheatmap)

# Assuming 'mME2' is available
# Create a new data frame with additional columns for land_use and season
# You might need to adjust the column names based on your actual data
metadata_df <- data.frame(
  sample = MEs0$sample,
  land_use = df_labels$Land_use,  # Assuming 'df_labels' contains land_use information
  season = df_labels$Season       # Assuming 'df_labels' contains season information
)

# Combine metadata with module data
combined_data <- cbind(metadata_df, mME)

heatmap_data <- combined_data[, c("sample", "name", "value", "land_use", "season")]

# Pivot the data to create a matrix for pheatmap
heatmap_matrix <- reshape2::dcast(heatmap_data, sample ~ name, value.var = "value")
pheatmap_metadata <- merge(metadata_df, heatmap_matrix)
# Create a matrix with values for pheatmap
pheatmap_data <- data.matrix(pheatmap_metadata[,-1:-3]) # Exclude non-numeric columns
rownames(pheatmap_data) <- pheatmap_metadata$sample
pheatmap_data <- t(pheatmap_data)
# Extract land_use and season columns for annotation
annotation_cols <- pheatmap_metadata[, c("sample","land_use", "season")]
annotation_df <- data.frame(pheatmap_metadata[,c("sample", "land_use", "season")])


# 
# install.packages("heatmaply")
library(heatmaply)
library(plotly)
heatmap <- heatmaply(
  pheatmap_data,
  Rowv = FALSE,
  Colv = FALSE,
  colors = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Module-Sample Relationships",
  xlab = "Samples",
  ylab = "Modules",
  k_col = 4,  # Number of clusters for columns (adjust as needed)
  layout_as_tree = TRUE,  # Show dendrogram as a tree
  add_heatmap_label = TRUE,
  col_side_colors = annotation_df[, c("land_use", "season")],  # Use col_side_colors for annotations on the columns
  fontsize_row = 4,  # Adjust the font size for rows
  fontsize_col = 4,
  file = "EN_Module_Sample_Relationship_heatmap_landuse_season.pdf"
)
# orca(heatmap, file = "heatmaply_heatmap.pdf")
pdf("EN_Module_Sample_Relationship_heatmap_landuse_season.pdf", width = 8, height = 6)
# print(heatmap)
# dev.off()



#Finding correlations between module-sample and physicochemical properties 
sample_modules <- read.delim("sample_modules.txt", header = TRUE, stringsAsFactors = FALSE)
file_path <- "Metadata_clean_AMRVFMG_project.csv"
metadata <-read.csv(file_path, header = TRUE, row.names = 1)
#figure out columns and their names to excludethem in the smaller dataframe

selected_columns <- c("Sample_Site", "sampleID", "date", "Season", "Land_use", "sample_type", "site", "lon", "lat", "Year", "Week", "yeardayofyear", "desc_", "site_yr")
pcp_dataframe <- metadata[selected_columns]
metadata_column_names <- colnames(pcp_dataframe)
                                  

pcp_df <- metadata[, !colnames(metadata) %in% metadata_column_names]


pcp_module <- merge(pcp_df, sample_modules, by.x = 0, by.y = "genus", all.x = TRUE)
# 
# install.packages("corrplot")
# library(corrplot)


# Convert the "cyl" variable to a factor
pcp_module$colors <- as.factor(pcp_module$colors)


# Filter rows with complete data for the 'colors' column
MEs0 <- moduleEigengenes(filtered_taxonomic_data_numeric_sample, mergedColors)$eigengenes
#reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)
pcp_module_filtered <- pcp_module[complete.cases(pcp_module$colors), ]
remove_columns <- c("Row.names", "site_type", "AMIA_AMN_F", "NITRITE_F", "NITRATE_F", "TOTKN_F", "TOTPHO_F", "DOC_F", "colors")
pcp_module_filtered <- pcp_module_filtered[, !colnames(pcp_module_filtered) %in% remove_columns]
remove_landuse <- grep("_\\d{4}$", names(pcp_module_filtered))
pcp_module_filtered_1 <- pcp_module_filtered[, -remove_landuse]

# gene_modules <- netwk$colors
# eigengenes = moduleEigengenes(filtered_taxonomic_data_numeric, netwk$colors)
# module_trait_cor <- cor(gene_modules, pcp_module_filtered)
module_trait_cor <- cor(MEs0, pcp_module_filtered_1)
p_values <- corPvalueStudent(module_trait_cor, nSamples = ncol(MEs0))



library(heatmaply)

# Create a heatmap using heatmaply
heatmaply(
  module_trait_cor,
  labRow = row.names(module_trait_cor),
  labCol = colnames(module_trait_cor),
  colors = colorRampPalette(c("blue", "white", "red"))(100),
  k_col = 4,  # Number of clusters for columns (adjust as needed)
  layout_as_tree = TRUE,  # Show dendrogram as a tree
  add_heatmap_label = TRUE,
  main = "Heatmap of MEcolor-Variable Correlations",
  xlab = "Variables",
  ylab = "MEcolors",
  Rowv = TRUE,
  Colv = TRUE,
  fontsize_row = 4,  # Adjust the font size for rows
  fontsize_col = 8,
  file = "EN_heatmap_MEcolor-variable corrlations.png", width = 800, height = 600
)

#filter correlations that have an absolute value of 0.3 of higher
module_trait_cor_0.3 <- module_trait_cor[abs(module_trait_cor) <= 0.3] <- 0
# make a dataframe of just the significant correlations
correlation_data <- data.frame(
  ME = rep(rownames(module_trait_cor), each = ncol(module_trait_cor)),
  Physicochemical_Parameter = rep(colnames(module_trait_cor), times = nrow(module_trait_cor)),
  Correlation = as.vector(module_trait_cor)
)
correlation_data <- as.data.frame(correlation_data)
# Filter based on the absolute value of correlation being greater than 0.3
correlation_data_filtered <- correlation_data[abs(correlation_data$Correlation) > 0.3, ]
correlation_data_filtered <- correlation_data_filtered[order(correlation_data_filtered$ME), ]
write.csv(correlation_data_filtered, file= "Correlation_significant_ME_PCP.csv", row.names = FALSE)



# Print the resulting matrix or data frame
print(correlation_matrix)


#Getting Pvalues for significant correlations

correlation_pvalue_data <- data.frame(
  ME = rep(rownames(p_values), each = ncol(p_values)),
  Physicochemical_Parameter = rep(colnames(p_values), times = nrow(p_values)),
  Pvalues = as.vector(p_values)
)
correlation_pvalue_data <- as.data.frame(correlation_pvalue_data)
correlation_pvalue_data_filtered <- correlation_pvalue_data[abs(correlation_pvalue_data$Pvalues) < 0.05, ]
correlation_pvalue_data_filtered <- correlation_pvalue_data_filtered[order(correlation_pvalue_data_filtered$ME), ]
write.csv(correlation_pvalue_data_filtered, file= "Correlation_significantpvalue_ME_PCP.csv", row.names = FALSE)

#####land_use and Season and MEs####
MEs0 <- moduleEigengenes(filtered_taxonomic_data_numeric_sample, mergedColors)$eigengenes
#reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)
filter_taxonomic_dataset <- t(as.data.frame(filtered_taxonomic_data))  # Ensure it's a data frame

# Check the resulting data frame
class(filter_taxonomic_dataset)
head(filter_taxonomic_dataset)  # Print the first few rows to verify the structure

# Select columns "Season" and "Land_use"
selected_columns <- c("Season", "Land_use")
metadata_dataframe <- filter_taxonomic_dataset[, selected_columns]

# module_LUS_cor <- cor(metadata_dataframe, MEs0)

# Combine the dataframes
combined_df <- cbind(metadata_dataframe, MEs0)
metadata <- combined_df[, c("Land_use", "Season")]
##################ANOVA and Tukey for all MES and Land_USE and Season#################################
# Create an empty list to store ANOVA results
anova_results_list <- list()

# Loop over each numeric column (ME) in MEs0
# Loop over each numeric column (ME) in MEs0
for (numeric_col in names(MEs0)) {
  # Perform ANOVA for each numeric column and both categorical variables
  anova_results_land_use <- aov(MEs0[[numeric_col]] ~ Land_use, data = combined_df)
  anova_results_season <- aov(MEs0[[numeric_col]] ~ Season, data = combined_df)
  
  # Perform Tukey's HSD for pairwise comparisons
  tukey_land_use <- TukeyHSD(anova_results_land_use)
  tukey_season <- TukeyHSD(anova_results_season)
  
  # Store the ME name, ANOVA, and pairwise comparison results in the list
  anova_results_list[[paste(numeric_col, "vs Land Use")]] <- list(
    ME = numeric_col, 
    ANOVA = summary(anova_results_land_use),
    Tukey = tukey_land_use
  )
  
  anova_results_list[[paste(numeric_col, "vs Season")]] <- list(
    ME = numeric_col, 
    ANOVA = summary(anova_results_season),
    Tukey = tukey_season
  )
}

# Print or further process the results in the list
for (result in anova_results_list) {
  cat("ME Name:", result$ME, "\n")
  cat("ANOVA Results:\n")
  print(result$ANOVA)
  cat("Tukey HSD Results:\n")
  print(result$Tukey)
}
#####Given overall correlation for MEs with Land_use and Season######################
metadata$Season <- as.numeric(as.factor(metadata$Season))
metadata$Land_use <-as.numeric(as.factor(metadata$Land_use))

module_LUS_cor <- cor(metadata, MEs0)

heatmaply(
  t(module_LUS_cor),
  labCol = row.names(module_LUS_cor),
  labRow = colnames(module_LUS_cor),
  colors = colorRampPalette(c("blue", "white", "red"))(100),
  k_col = 4,  # Number of clusters for columns (adjust as needed)
  layout_as_tree = TRUE,  # Show dendrogram as a tree
  add_heatmap_label = TRUE,
  main = "Heatmap of MEcolor-Variable Correlations",
  xlab = "Variables",
  ylab = "MEcolors",
  Rowv = FALSE,
  Colv = FALSE,
  fontsize_row = 4,  # Adjust the font size for rows
  fontsize_col = 8,
  file = "EN_heatmap_MEcolor-variable correlations_landuseseason.png", width = 800, height = 600
)

#####################correlation between ME and levels of land_use and season
combined_df$Season <- as.numeric(as.factor(combined_df$Season))
combined_df$Land_use <- as.numeric(as.factor(combined_df$Land_use))

numeric_columns <- names(combined_df)[sapply(combined_df, is.numeric) & !(names(combined_df) %in% c('Land_use', 'Season'))]

correlation_df <- data.frame(Land_use = character(),
                             Season = character(),
                             stringsAsFactors = FALSE)

# Loop over numeric columns
for (numeric_col in numeric_columns) {
  # Calculate correlations for each level of 'Land_use'
  correlations_land_use <- by(combined_df[[numeric_col]], combined_df$Land_use, cor)
  
  # Calculate correlations for each level of 'Season'
  correlations_season <- by(combined_df[[numeric_col]], combined_df$Season, cor)
  
  # Extract correlation values and add to the dataframe
  for (land_use_level in names(correlations_land_use)) {
    correlation_df <- rbind(correlation_df,
                            data.frame(Land_use = land_use_level,
                                       Season = "All",
                                       ME = numeric_col,
                                       Correlation = correlations_land_use[[land_use_level]]))
  }
  
  for (season_level in names(correlations_season)) {
    correlation_df <- rbind(correlation_df,
                            data.frame(Land_use = "All",
                                       Season = season_level,
                                       ME = numeric_col,
                                       Correlation = correlations_season[[season_level]]))
  }
}

# Print or further process the correlation dataframe
print(correlation_df)

df_combined <- correlation_df %>%
  mutate(Combined = paste(Land_use, Season, sep = "_"))


df_wide <- df_combined %>%
  pivot_wider(names_from = ME, values_from = Correlation)

df_wide_combinedname <- df_wide[, !(colnames(df_wide) %in% c("Land_use", "Season"))]
df_noname <- df_wide_combinedname[, !(colnames(df_wide_combinedname) %in% c("Combined"))]


heatmaply(
  t(df_noname),
  labCol = c("Mixed", "Agri-ditch", "Forested", "Spring", "Summer", "Autumn"),
  #labRow = colnames(df_wide_combinedname),
  colors = colorRampPalette(c("blue", "white", "red"))(100),
  k_col = 4,  # Number of clusters for columns (adjust as needed)
  layout_as_tree = TRUE,  # Show dendrogram as a tree
  add_heatmap_label = TRUE,
  main = "Heatmap of MEcolor-Variable Correlations",
  xlab = "Variables",
  ylab = "MEcolors",
  Rowv = FALSE,
  Colv = FALSE,
  fontsize_row = 4,  # Adjust the font size for rows
  fontsize_col = 8,
  file = "EN_heatmap_MEcolor-variable correlations_landuseseason_separateclasses.png", width = 800, height = 600 
)
##################General Linear Mixed Effects Models to check the correlation relationship between modules/members with the selected env. variables#############################
#picked modules that have the strong (>0.5) correlations: MEbrown4, MEdarkorange2, MEdarkslateblue, MElightslateblue; MElightsteelblue1, MElightyellow, MEmagenta, MEmidnightblue, MEorange, MEplum4, MEsalmon1, MEtan, MEturqoise

###making the dataframe for the analysis#############
#ME values + metadata
file_path <- "Metadata_clean_AMRVFMG_project.csv"
metadata_glme <-read.csv(file_path, header = TRUE, row.names = 1)
#figure out columns and their names to excludethem in the smaller dataframe

selected_columns <- c("sampleID", "sample_type", "site", "lon", "lat", "Year", "Week", "yeardayofyear", "desc_", "site_yr")

metadata_filtered <- metadata_glme[,!(colnames(metadata_glme) %in% selected_columns)]

MEs0$Sample_ID = row.names(t(filtered_taxonomic_data))
metadata_filtered$Sample_ID <- rownames(metadata_filtered)

merged_df_glme <- merge(MEs0, metadata_filtered, by = 'Sample_ID')


#create blocks
merged_df_glme$Block <- ifelse(merged_df_glme$Sample_Site %in% c("SN018", "SN019"), "Block1819",
                                       ifelse(merged_df_glme$Sample_Site %in% c("SN020", "SN021"), "Block2021", paste0("Block_", merged_df_glme$Sample_Site)))
#convert date to a numeric
# Assuming your "date" variable is in Date format, convert it to numeric (days since a reference date)
merged_df_glme$date <- as.Date(merged_df_glme$date, format = "%m/%d/%y")
merged_df_glme$date_numeric <- as.numeric(merged_df_glme$date - min(merged_df_glme$date))

################GLME, print to file
#Get packages for glme
#install.packages("lme4")
#install,packages("car")
library(lme4)
library(car)

#for MEbrown

file_path <- "glme_contdate_MEbrown.txt"
response_variables <- c("strahler", "PH", "CONDUCTIVITY_MSC", "DISS_OXYGEN_P", "DISS_OXYGEN_MGL", "ORP_MV", "rain_mm_2d", "rain_mm_7d", "avg_temp_C_1d", "avg_temp_C_2d", "avg_temp_C_3d", "RU_DISM3S", "BE_DISM3S")


capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with lme4
    formula_str <- paste("MEbrown ~", response_variable, "+ (1 | Block) + (1 | date)", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      # Check for missing values in MEbrown4
      if (any(is.na(merged_df_glme$MEbrown))) {
        cat("Error: Missing values in predictor variable (MEbrown4) for", response_variable, "\n")
      } else {
        model <- glmer(formula,
                       data = merged_df_glme,
                       family = gaussian(link = "identity"))
        
        # Print results for the current response variable
        cat("Response Variable:", response_variable, "\n")
        print(summary(model))
        
        cat("ANOVA:\n")
        anova_result <- Anova(model, type = "II")
        print(anova_result)
        
        # Extract p-value from the ANOVA result
        p_value <- anova_result$`Pr(>Chisq)`[2]
        cat("P-Value:", p_value, "\n\n")
      }
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
  
  # Close the file connection
  sink()
}, file = file_path, append = TRUE)

# Check the file contents
cat(readLines(file_path), sep = "\n")



#For MEdarkmagenta

response_variables <- c("TEMP_C",
                        "NITRATE",
                        "TOTKN",
                        "TOTKN_MDL",
                        "TOTPHO",
                        "TOTPHO_MDL")

# Open a file connection for writing
file_path <- "glme_contdate_MEdarkmagenta.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with lme4
    formula_str <- paste("MEdarkmagenta ~", response_variable, "+ (1 | Block) + (1 | date)", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      # Check for missing values in ME
      if (any(is.na(merged_df_glme$MEdarkmagenta))) {
        cat("Error: Missing values in predictor variable  for", response_variable, "\n")
      } else {
        model <- glmer(formula,
                       data = merged_df_glme,
                       family = gaussian(link = "identity"))
        
        # Print results for the current response variable
        cat("Response Variable:", response_variable, "\n")
        print(summary(model))
        
        cat("ANOVA:\n")
        anova_result <- Anova(model, type = "II")
        print(anova_result)
        
        # Extract p-value from the ANOVA result
        p_value <- anova_result$`Pr(>Chisq)`[2]
        cat("P-Value:", p_value, "\n\n")
      }
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
  
  # Close the file connection
  sink()
}, file = file_path, append = TRUE)

# Check the file contents
cat(readLines(file_path), sep = "\n")


#For MEdarkorange2

response_variables <- c("DISS_OXYGEN_MGL",
                        "TURBIDITY_NTU",
                        "AMIA_AMN",
                        "AMIA_AMN_MDL",
                        "NITRITE",
                        "NITRITE_MDL",
                        "avg_temp_C_2d",
                        "avg_temp_C_7d",
                        "min_temp_C_1d",
                        "max_temp_C_1d",
                        "RU_DISM3S")

# Open a file connection for writing
file_path <- "glme_contdate_MEdarkorange2.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with lme4
    formula_str <- paste("MEdarkorange2 ~", response_variable, "+ (1 | Block) + (1 | date)", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      # Check for missing values in ME
      if (any(is.na(merged_df_glme$MEdarkorange2))) {
        cat("Error: Missing values in predictor variable for", response_variable, "\n")
      } else {
        model <- glmer(formula,
                       data = merged_df_glme,
                       family = gaussian(link = "identity"))
        
        # Print results for the current response variable
        cat("Response Variable:", response_variable, "\n")
        print(summary(model))
        
        cat("ANOVA:\n")
        anova_result <- Anova(model, type = "II")
        print(anova_result)
        
        # Extract p-value from the ANOVA result
        p_value <- anova_result$`Pr(>Chisq)`[2]
        cat("P-Value:", p_value, "\n\n")
      }
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
  
  # Close the file connection
  sink()
}, file = file_path, append = TRUE)

# Check the file contents
cat(readLines(file_path), sep = "\n")



#For MEdarkturquoise

response_variables <- c("TOC_MDL",
                        "rain_mm_7d",
                        "min_temp_C_1d",
                        "RU_DISM3S",
                        "BE_DISM3S",
                        "daily_avg_solar_radiation_wm2")

# Open a file connection for writing
file_path <- "glme_contdate_MEdarkturquoise.txt"
capture.output({
  sink(file_path, append = FALSE)

  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with lme4
    formula_str <- paste("MEdarkturquoise ~", response_variable, "+ (1 | Block) + (1 | date)", sep = "")
    formula <- as.formula(formula_str)

    tryCatch({
      # Check for missing values in ME
      if (any(is.na(merged_df_glme$MEdarkturquoise))) {2
        cat("Error: Missing values in predictor variable for", response_variable, "\n")
      } else {
        model <- glmer(formula,
                       data = merged_df_glme,
                       family = gaussian(link = "identity"))

        # Print results for the current response variable
        cat("Response Variable:", response_variable, "\n")
        print(summary(model))

        cat("ANOVA:\n")
        anova_result <- Anova(model, type = "II")
        print(anova_result)

        # Extract p-value from the ANOVA result
        p_value <- anova_result$`Pr(>Chisq)`[2]
        cat("P-Value:", p_value, "\n\n")
      }
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }

  # Close the file connection
  sink()
}, file = file_path, append = TRUE)

# Check the file contents
cat(readLines(file_path), sep = "\n")


#For MEgreen

response_variables <- c("ORP_MV",
                        "AMIA_AMN_MDL",
                        "NITRITE",
                        "NITRITE_MDL",
                        "NITRATE",
                        "TOTPHO_MDL",
                        "DOC",
                        "TOC_MDL")

# Open a file connection for writing
file_path <- "glme_contdate_green.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with lme4
    formula_str <- paste("MEgreen ~", response_variable, "+ (1 | Block) + (1 | date)", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      # Check for missing values in ME
      if (any(is.na(merged_df_glme$MEgreen))) {
        cat("Error: Missing values in predictor variable for", response_variable, "\n")
      } else {
        model <- glmer(formula,
                       data = merged_df_glme,
                       family = gaussian(link = "identity"))
        
        # Print results for the current response variable
        cat("Response Variable:", response_variable, "\n")
        print(summary(model))
        
        cat("ANOVA:\n")
        anova_result <- Anova(model, type = "II")
        print(anova_result)
        
        # Extract p-value from the ANOVA result
        p_value <- anova_result$`Pr(>Chisq)`[2]
        cat("P-Value:", p_value, "\n\n")
      }
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
  
  # Close the file connection
  sink()
}, file = file_path, append = TRUE)

# Check the file contents
cat(readLines(file_path), sep = "\n")

#For MEivory

response_variables <- c("TOTPHO_MDL",
                        "TOC_MDL",
                        "avg_temp_C_1d",
                        "avg_temp_C_3d",
                        "avg_temp_C_5d",
                        "avg_temp_C_7d",
                        "min_temp_C_1d",
                        "max_temp_C_1d")

# Open a file connection for writing
file_path <- "glme_contdate_MEivory.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with lme4
    formula_str <- paste("MEivory ~", response_variable, "+ (1 | Block) + (1 | date)", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      # Check for missing values in ME
      if (any(is.na(merged_df_glme$MEivory))) {
        cat("Error: Missing values in predictor variable for", response_variable, "\n")
      } else {
        model <- glmer(formula,
                       data = merged_df_glme,
                       family = gaussian(link = "identity"))
        
        # Print results for the current response variable
        cat("Response Variable:", response_variable, "\n")
        print(summary(model))
        
        cat("ANOVA:\n")
        anova_result <- Anova(model, type = "II")
        print(anova_result)
        
        # Extract p-value from the ANOVA result
        p_value <- anova_result$`Pr(>Chisq)`[1]
        cat("P-Value:", p_value, "\n\n")
      }
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
  
  # Close the file connection
  sink()
}, file = file_path, append = TRUE)

# Check the file contents
cat(readLines(file_path), sep = "\n")

#for MElightcyan
response_variables <- c("ORP_MV",
                        "DOC_MDL",
                        "TOC",
                        "TOC_F",
                        "TOC_MDL")

# Open a file connection for writing
file_path <- "glme_contdate_MElightcyan.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with lme4
    formula_str <- paste("MElightcyan ~", response_variable, "+ (1 | Block) + (1 | date)", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      # Check for missing values in ME
      if (any(is.na(merged_df_glme$MElightcyan))) {
        cat("Error: Missing values in predictor variable for", response_variable, "\n")
      } else {
        model <- glmer(formula,
                       data = merged_df_glme,
                       family = gaussian(link = "identity"))
        
        # Print results for the current response variable
        cat("Response Variable:", response_variable, "\n")
        print(summary(model))
        
        cat("ANOVA:\n")
        anova_result <- Anova(model, type = "II")
        print(anova_result)
        
        # Extract p-value from the ANOVA result
        p_value <- anova_result$`Pr(>Chisq)`[1]
        cat("P-Value:", p_value, "\n\n")
      }
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
  
  # Close the file connection
  sink()
}, file = file_path, append = TRUE)

# Check the file contents
cat(readLines(file_path), sep = "\n")

#for MElightcyan1
response_variables <- c("DISS_OXYGEN_P",
                        "NITRITE",
                        "NITRATE",
                        "NITRATE_MDL",
                        "TOTKN",
                        "TOTKN_MDL",
                        "TOTPHO",
                        "max_temp_C_1d",
                        "daily_avg_solar_radiation_wm2")

# Open a file connection for writing
file_path <- "glme_contdate_MElightcyan1.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with lme4
    formula_str <- paste("MElightcyan1 ~", response_variable, "+ (1 | Block) + (1 | date)", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      # Check for missing values in ME
      if (any(is.na(merged_df_glme$MElightcyan1))) {
        cat("Error: Missing values in predictor variable for", response_variable, "\n")
      } else {
        model <- glmer(formula,
                       data = merged_df_glme,
                       family = gaussian(link = "identity"))
        
        # Print results for the current response variable
        cat("Response Variable:", response_variable, "\n")
        print(summary(model))
        
        cat("ANOVA:\n")
        anova_result <- Anova(model, type = "II")
        print(anova_result)
        
        # Extract p-value from the ANOVA result
        p_value <- anova_result$`Pr(>Chisq)`[1]
        cat("P-Value:", p_value, "\n\n")
      }
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
  
  # Close the file connection
  sink()
}, file = file_path, append = TRUE)

# Check the file contents
cat(readLines(file_path), sep = "\n")

#for MElightsteelblue
response_variables <- c("TEMP_C",
                        "TURBIDITY_NTU",
                        "AMIA_AMN_MDL",
                        "NITRITE",
                        "NITRITE_MDL",
                        "NITRATE",
                        "NITRATE_MDL",
                        "avg_temp_C_5d",
                        "max_temp_C_1d",
                        "RU_DISM3S",
                        "BE_DISM3S",
                        "daily_avg_solar_radiation_wm2")

# Open a file connection for writing
file_path <- "glme_contdate_MElightsteelblue.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with lme4
    formula_str <- paste("MElightsteelblue ~", response_variable, "+ (1 | Block) + (1 | date)", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      # Check for missing values in ME
      if (any(is.na(merged_df_glme$MElightsteelblue))) {
        cat("Error: Missing values in predictor variable for", response_variable, "\n")
      } else {
        model <- glmer(formula,
                       data = merged_df_glme,
                       family = gaussian(link = "identity"))
        
        # Print results for the current response variable
        cat("Response Variable:", response_variable, "\n")
        print(summary(model))
        
        cat("ANOVA:\n")
        anova_result <- Anova(model, type = "II")
        print(anova_result)
        
        # Extract p-value from the ANOVA result
        p_value <- anova_result$`Pr(>Chisq)`[1]
        cat("P-Value:", p_value, "\n\n")
      }
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
  
  # Close the file connection
  sink()
}, file = file_path, append = TRUE)

# Check the file contents
cat(readLines(file_path), sep = "\n")

#for MElightyellow
response_variables <- c("strahler",
                        "TEMP_C",
                        "PH",
                        "TURBIDITY_NTU",
                        "AMIA_AMN",
                        "NITRATE")

# Open a file connection for writing
file_path <- "glme_contdate_MElightyellow.txt"
capture.output({
  sink(file_path, append = FALSE)

  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with lme4
    formula_str <- paste("MElightyellow ~", response_variable, "+ (1 | Block) + (1 | date)", sep = "")
    formula <- as.formula(formula_str)

    tryCatch({
      # Check for missing values in ME
      if (any(is.na(merged_df_glme$MElightyellow))) {
        cat("Error: Missing values in predictor variable for", response_variable, "\n")
      } else {
        model <- glmer(formula,
                       data = merged_df_glme,
                       family = gaussian(link = "identity"))

        # Print results for the current response variable
        cat("Response Variable:", response_variable, "\n")
        print(summary(model))

        cat("ANOVA:\n")
        anova_result <- Anova(model, type = "II")
        print(anova_result)

        # Extract p-value from the ANOVA result
        p_value <- anova_result$`Pr(>Chisq)`[1]
        cat("P-Value:", p_value, "\n\n")
      }
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }

  # Close the file connection
  sink()
}, file = file_path, append = TRUE)

# Check the file contents
cat(readLines(file_path), sep = "\n")

#for MEmediumpurple2
response_variables <- c("NITRATE_MDL",
                        "rain_mm_5d",
                        "rain_mm_7d",
                        "avg_temp_C_1d",
                        "avg_temp_C_2d",
                        "avg_temp_C_3d")

# Open a file connection for writing
file_path <- "glme_contdate_MEmediumpurple2.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with lme4
    formula_str <- paste("MEmediumpurple2 ~", response_variable, "+ (1 | Block) + (1 | date)", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      # Check for missing values in MEdarkorange2
      if (any(is.na(merged_df_glme$MEmediumpurple2))) {
        cat("Error: Missing values in predictor variable (MEdarkorange2) for", response_variable, "\n")
      } else {
        model <- glmer(formula,
                       data = merged_df_glme,
                       family = gaussian(link = "identity"))
        
        # Print results for the current response variable
        cat("Response Variable:", response_variable, "\n")
        print(summary(model))
        
        cat("ANOVA:\n")
        anova_result <- Anova(model, type = "II")
        print(anova_result)
        
        # Extract p-value from the ANOVA result
        p_value <- anova_result$`Pr(>Chisq)`[1]
        cat("P-Value:", p_value, "\n\n")
      }
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
  
  # Close the file connection
  sink()
}, file = file_path, append = TRUE)

# Check the file contents
cat(readLines(file_path), sep = "\n")

#for MEorangered3
response_variables <- c("AMIA_AMN_MDL",
                        "TOTPHO_MDL",
                        "DOC",
                        "DOC_MDL",
                        "TOC",
                        "TOC_F",
                        "avg_temp_C_2d")

# Open a file connection for writing
file_path <- "glme_contdate_MEorangered3.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with lme4
    formula_str <- paste("MEorangered3 ~", response_variable, "+ (1 | Block) + (1 | date)", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      # Check for missing values in MEdarkorange2
      if (any(is.na(merged_df_glme$MEorangered3))) {
        cat("Error: Missing values in predictor variable (MEdarkorange2) for", response_variable, "\n")
      } else {
        model <- glmer(formula,
                       data = merged_df_glme,
                       family = gaussian(link = "identity"))
        
        # Print results for the current response variable
        cat("Response Variable:", response_variable, "\n")
        print(summary(model))
        
        cat("ANOVA:\n")
        anova_result <- Anova(model, type = "II")
        print(anova_result)
        
        # Extract p-value from the ANOVA result
        p_value <- anova_result$`Pr(>Chisq)`[1]
        cat("P-Value:", p_value, "\n\n")
      }
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
  
  # Close the file connection
  sink()
}, file = file_path, append = TRUE)

# Check the file contents
cat(readLines(file_path), sep = "\n")


#for MEplum1
response_variables <- c("NITRITE_MDL",
                        "NITRATE_MDL",
                        "TOTKN",
                        "TOTKN_MDL",
                        "TOTPHO",
                        "rain_mm_5d")

# Open a file connection for writing
file_path <- "glme_contdate_MEplum1.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with lme4
    formula_str <- paste("MEplum1 ~", response_variable, "+ (1 | Block) + (1 | date)", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      # Check for missing values in ME
      if (any(is.na(merged_df_glme$MEplum1))) {
        cat("Error: Missing values in predictor variable for", response_variable, "\n")
      } else {
        model <- glmer(formula,
                       data = merged_df_glme,
                       family = gaussian(link = "identity"))
        
        # Print results for the current response variable
        cat("Response Variable:", response_variable, "\n")
        print(summary(model))
        
        cat("ANOVA:\n")
        anova_result <- Anova(model, type = "II")
        print(anova_result)
        
        # Extract p-value from the ANOVA result
        p_value <- anova_result$`Pr(>Chisq)`[1]
        cat("P-Value:", p_value, "\n\n")
      }
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
  
  # Close the file connection
  sink()
}, file = file_path, append = TRUE)

# Check the file contents
cat(readLines(file_path), sep = "\n")


#for MEsaddlebrown
response_variables <- c("TOTKN",
                        "TOTPHO_MDL",
                        "DOC",
                        "DOC_MDL",
                        "TOC",
                        "rain_mm_3d",
                        "rain_mm_5d",
                        "avg_temp_C_3d")

# Open a file connection for writing
file_path <- "glme_contdate_MEsaddlebrown.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with lme4
    formula_str <- paste("MEsaddlebrown ~", response_variable, "+ (1 | Block) + (1 | date)", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      # Check for missing values in ME
      if (any(is.na(merged_df_glme$MEsaddlebrown))) {
        cat("Error: Missing values in predictor variable for", response_variable, "\n")
      } else {
        model <- glmer(formula,
                       data = merged_df_glme,
                       family = gaussian(link = "identity"))
        
        # Print results for the current response variable
        cat("Response Variable:", response_variable, "\n")
        print(summary(model))
        
        cat("ANOVA:\n")
        anova_result <- Anova(model, type = "II")
        print(anova_result)
        
        # Extract p-value from the ANOVA result
        p_value <- anova_result$`Pr(>Chisq)`[1]
        cat("P-Value:", p_value, "\n\n")
      }
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
  
  # Close the file connection
  sink()
}, file = file_path, append = TRUE)

# Check the file contents
cat(readLines(file_path), sep = "\n")



#for MEthistle
response_variables <- c("AMIA_AMN_MDL",
                        "TOTKN",
                        "TOC",
                        "TOC_MDL",
                        "rain_mm_1d",
                        "rain_mm_2d",
                        "rain_mm_3d",
                        "rain_mm_5d")

# Open a file connection for writing
file_path <- "glme_contdate_MEthistle.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with lme4
    formula_str <- paste("MEthistle ~", response_variable, "+ (1 | Block) + (1 | date)", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      # Check for missing values in ME
      if (any(is.na(merged_df_glme$MEthistle))) {
        cat("Error: Missing values in predictor variable for", response_variable, "\n")
      } else {
        model <- glmer(formula,
                       data = merged_df_glme,
                       family = gaussian(link = "identity"))
        
        # Print results for the current response variable
        cat("Response Variable:", response_variable, "\n")
        print(summary(model))
        
        cat("ANOVA:\n")
        anova_result <- Anova(model, type = "II")
        print(anova_result)
        
        # Extract p-value from the ANOVA result
        p_value <- anova_result$`Pr(>Chisq)`[1]
        cat("P-Value:", p_value, "\n\n")
      }
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
  
  # Close the file connection
  sink()
}, file = file_path, append = TRUE)

# Check the file contents
cat(readLines(file_path), sep = "\n")


#for MEthistle1
response_variables <- c("strahler",
                        "TEMP_C",
                        "DOC_MDL",
                        "TOC_MDL",
                        "rain_mm_1d",
                        "rain_mm_2d",
                        "rain_mm_3d",
                        "avg_temp_C_3d",
                        "avg_temp_C_5d")

# Open a file connection for writing
file_path <- "glme_contdate_MEthistle1.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with lme4
    formula_str <- paste("MEthistle1 ~", response_variable, "+ (1 | Block) + (1 | date)", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      # Check for missing values in ME
      if (any(is.na(merged_df_glme$MEthistle1))) {
        cat("Error: Missing values in predictor variable for", response_variable, "\n")
      } else {
        model <- glmer(formula,
                       data = merged_df_glme,
                       family = gaussian(link = "identity"))
        
        # Print results for the current response variable
        cat("Response Variable:", response_variable, "\n")
        print(summary(model))
        
        cat("ANOVA:\n")
        anova_result <- Anova(model, type = "II")
        print(anova_result)
        
        # Extract p-value from the ANOVA result
        p_value <- anova_result$`Pr(>Chisq)`[1]
        cat("P-Value:", p_value, "\n\n")
      }
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
  
  # Close the file connection
  sink()
}, file = file_path, append = TRUE)

# Check the file contents
cat(readLines(file_path), sep = "\n")



