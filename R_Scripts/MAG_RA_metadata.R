# Relative abundances of MAGs in metagenomic samples analysis
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1

library(tidyr)
library(dplyr)
#set your working directory were everything will be stored and saved
setwd("C:/Users/Sandersonh/OneDrive - AGR-AGR/Desktop/Wen Chen/Project_1_landuse_AMR_virulence/Statistics_Analysis")
MAG_counts <- read.csv("MAG_read_mapping_counts_table.csv", header = TRUE, stringsAsFactors = FALSE)
metadata <- read.csv("Metadata_clean_AMRVFMG_project.csv", header = TRUE, stringsAsFactors = FALSE)

# 
# transformed_MAG_counts <- MAG_counts %>%
#   pivot_longer(cols = starts_with("metabat_bins"), 
#                names_to = "MAG", 
#                values_to = "Read_Count")
#################Normalize MAG read counts
column_sums <- colSums(MAG_counts[, -1])  # Exclude the first column (Sample)

# Divide each value in the dataframe by its column sum
normalized_MAG_counts <- MAG_counts %>%
  mutate(across(starts_with("metabat_bins"), ~ . / column_sums))


#################clr transformation
library(compositions)
sample_names <- normalized_MAG_counts$SAMPLE_ID
MAG_data <- normalized_MAG_counts[, -which(colnames(normalized_MAG_counts) == "SAMPLE_ID")]
#CLR tranformation for Genus taxonomic data to generate a dataframe
# Perform CLR transformation
clr_data_MAG <- clr(MAG_data)

# Convert CLR-transformed data to a data frame
clr_data_MAG_df <- as.data.frame(clr_data_MAG)

# Add sample names as a column
clr_data_MAG_df$SAMPLE_ID <- sample_names

# Reorder columns (optional, adjust as needed)
clr_data_MAG_df <- clr_data_MAG_df[, c("SAMPLE_ID", colnames(clr_data_MAG_df)[1:(ncol(clr_data_MAG_df)-1)])]

# Display the CLR-transformed data as a data frame
print(clr_data_MAG_df)

#create a dataframe with Land_use and Season with the MAG counts
selected_columns <- c("Sample_ID_Genome", "Land_use", "Season", "Sample_Site", "date")
metadata_dataframe <- metadata[selected_columns]

metadata_dataframe <- metadata_dataframe %>% rename(SAMPLE_ID = Sample_ID_Genome)

combined_df <- merge(metadata_dataframe, clr_data_MAG_df, by = "SAMPLE_ID", all = TRUE)

#creating BLOCKs
combined_df$Block <- ifelse(combined_df$Sample_Site %in% c("SN018", "SN019"), "Block1819",
                                       ifelse(combined_df$Sample_Site %in% c("SN020", "SN021"), "Block2021", paste0("Block_", combined_df$Sample_Site)))
#creating numerical dates
combined_df$date <- as.Date(combined_df$date, format = "%m/%d/%y")
combined_df$date_numeric <- as.numeric(combined_df$date - min(combined_df$date))

#response variable for each MAG column
combined_df <- na.omit(combined_df)
response_variables <- grep("^metabat_bins", names(combined_df), value = TRUE)
combined_df[response_variables] <- lapply(combined_df[response_variables], as.numeric)

#calculate LME for Land_Use and Season and MAG counts
library(nlme)
library(emmeans)
# Open a file connection for writing
file_path <- "MAG_lme_contdate_landuse_results.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ Land_use", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = combined_df)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      print(summary(model))
      
      # ANOVA
      cat("ANOVA:\n")
      print(anova(model))
      
      # Pairwise comparisons
      cat("Pairwise Comparison Results:\n")
      pwc <- emmeans(model, specs = "Land_use", data = combined_df)
      pairwise_results <- pairs(pwc, adjust = "BH")
      print(pairwise_results)
      
      cat("\n")
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


#calculate LME forSeason and MAG counts
library(nlme)
library(emmeans)
# Open a file connection for writing
file_path <- "MAG_lme_contdate_season_results.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ Season", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = combined_df)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      print(summary(model))
      
      # ANOVA
      cat("ANOVA:\n")
      print(anova(model))
      
      # Pairwise comparisons
      cat("Pairwise Comparison Results:\n")
      pwc <- emmeans(model, specs = "Season", data = combined_df)
      pairwise_results <- pairs(pwc, adjust = "BH")
      print(pairwise_results)
      
      cat("\n")
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

#LME MAGS and Interaction Effects
library(broom)
file_path2 <-"lme_contdate_seasonlanduse_results_pvalues.txt"

for (response_variable in response_variables) {
  # Fit the linear mixed-effects model with lme, including the interaction term
  formula_str <- paste(response_variable, "~ Land_use * Season", sep = "")
  formula <- as.formula(formula_str)
  
  tryCatch({
    model <- lme(fixed = formula,
                 random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                 data = combined_df)
    
    # Print results for the current response variable
    cat("Response Variable:", response_variable, "\n")
    
    # ANOVA
    cat("ANOVA:\n")
    anova_result <- anova(model)
    df_anova_result <- data.frame(anova_result)
    # Display ANOVA results
    print(df_anova_result)
    
    # Extract specific values from ANOVA results
    land_use_p_value <- df_anova_result["Land_use", "p.value"]
    season_p_value <- df_anova_result["Season", "p.value"]
    interaction_p_value <- df_anova_result["Land_use:Season", "p.value"]
    
    cat("Land_Use p-value: ", land_use_p_value, "\n")
    cat("Season p-value: ", season_p_value, "\n")
    cat("Interaction (Land_Use:Season) p-value: ", interaction_p_value, "\n")
    
    # Perform p-value correction
    land_use_p_value_corrected <- p.adjust(land_use_p_value, method = "bonferroni")
    season_p_value_corrected <- p.adjust(season_p_value, method = "bonferroni")
    interaction_p_value_corrected <- p.adjust(interaction_p_value, method = "bonferroni") 
    
    # Save ANOVA results to data frame
    anova_df <- rbind(anova_df, data.frame(Response_Variable = response_variable,
                                           Land_Use_PValue = land_use_p_value_corrected,
                                           Season_PValue = season_p_value_corrected,
                                           Interaction_PValue = interaction_p_value_corrected,
                                           stringsAsFactors = FALSE))
    
    
    cat("Land_Use p-value (corrected): ", land_use_p_value_corrected, "\n")
    cat("Season p-value (corrected): ", season_p_value_corrected, "\n")
    cat("Interaction (Land_Use:Season) p-value (corrected): ", interaction_p_value_corrected, "\n")
    
    # Pairwise comparisons
    cat("Pairwise Comparison Results:\n")
    pwc <- emmeans(model, specs = c("Land_use", "Season"))
    pairwise_results <- tidy(pwc, conf.int = TRUE)
    print(pairwise_results)
    
    # Save pairwise results to data frame
    pairwise_df <- rbind(pairwise_df, data.frame(Response_Variable = response_variable,
                                                 Comparison1 = pairwise_results[,1],
                                                 Comparison2 = pairwise_results[,2],
                                                 Estimate = pairwise_results$estimate,
                                                 SE = pairwise_results$std.error,
                                                 df = pairwise_results$df,
                                                 conf.low = pairwise_results$conf.low,
                                                 conf.high = pairwise_results$conf.high,
                                                 statistic = pairwise_results$statistic,
                                                 PValue = p.adjust(pairwise_results$p.value, method = "bonferroni"),
                                                 stringsAsFactors = FALSE))
    
    cat("\n")
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

# Print data frames
print(anova_df)
print(pairwise_df)

write.csv(anova_df, file = "MAG_lme_ANOVA_pvalues_seasonlanduse.csv", row.names = FALSE)
write.csv(pairwise_df, file = "MAG_lme_PWC_pvalues_seasonlanduse.csv", row.names = FALSE)
####Combine the two tables
anova_S <- read.csv("MAG_lme_ANOVA_pvalues_seasonlanduse.csv", header = TRUE, stringsAsFactors = FALSE) 
pairwise_S <- read.csv("MAG_lme_PWC_pvalues_seasonlanduse.csv", header = TRUE, stringsAsFactors = FALSE) 
merged_table_pvalues_S <- merge(anova_S, pairwise_S, by.x = 'Response_Variable', by.y = 'Response_Variable')
write.csv(merged_table_pvalues_S, file = "MAG_lme_ALL_pvalues_seasonlanduse.csv", row.names = FALSE)


################Making figures for MAG counts Land_use
library(ggplot2)
library(ggpubr)
library(gridExtra)

metadata$Land_use <- factor(metadata$Land_use, levels = c("Agri_ditch", "Mixed", "Forested"))
land_use_colors <- c("Agri_ditch" = "black", "Mixed" = "green", "Forested" = "red")
plot_list_LU <- list()

#MAG 14
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_14_LU <- c("0.1799", "0.1553", "0.1799")


# Create boxplot
plot_list_LU[["14_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.14, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.14, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_14_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.14, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_14_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.14, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_14_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["14_LU"]])

#MAG 17

pairwise_p_values_17_LU <- c('0.9686','0.9686','0.9686')


# Create boxplot
plot_list_LU[["17_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.17, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.17, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_17_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.17, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_17_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.17, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_17_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["17_LU"]])


#mag 42

pairwise_p_values_42_LU <- c('0.0164','0.5859','0.0164')


# Create boxplot
plot_list_LU[["42_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.42, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.42, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_42_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.42, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_42_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.42, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_42_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["42_LU"]])

#mag 47
pairwise_p_values_47_LU <- c('0.1502','0.5452','0.1502')


# Create boxplot
plot_list_LU[["47_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.47, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.47, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_47_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.47, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_47_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.47, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_47_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["47_LU"]])


#mag 166

pairwise_p_values_166_LU <- c('0.7473','0.2423','0.1655')


# Create boxplot
plot_list_LU[["166_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.166, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.166, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_166_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.166, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_166_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.166, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_166_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["166_LU"]])

#mag 209
pairwise_p_values_209_LU <- c('0.6295','0.6295','0.9797')


# Create boxplot
plot_list_LU[["209_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.209, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.209, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_209_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.209, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_209_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.209, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_209_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["209_LU"]])


#mag 240
pairwise_p_values_240_LU <- c('0.8691','0.3157','0.2574')


# Create boxplot
plot_list_LU[["240_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.240, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.240, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_240_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.240, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_240_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.240, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_240_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["240_LU"]])

#mag 283
pairwise_p_values_283_LU <- c('0.6591','0.6591','0.9412')


# Create boxplot
plot_list_LU[["283_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.283, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.283, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_283_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.283, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_283_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.283, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_283_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["283_LU"]])

#mag 369
pairwise_p_values_369_LU <- c('0.1903','0.6243','0.1903')


# Create boxplot
plot_list_LU[["369_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.369, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.369, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_369_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.369, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_369_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.369, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_369_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["369_LU"]])


#mag 530
pairwise_p_values_530_LU <- c('0.2635','0.1166','0.0329')


# Create boxplot
plot_list_LU[["530_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.530, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.530, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_530_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.530, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_530_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.530, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_530_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["530_LU"]])


#mag 616
pairwise_p_values_616_LU <- c('0.6465','0.3799','0.2259')


# Create boxplot
plot_list_LU[["616_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.616, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.616, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_616_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.616, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_616_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.616, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_616_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["616_LU"]])


#mag 647

pairwise_p_values_647_LU <- c('0.0717','0.0717','0.8309')


# Create boxplot
plot_list_LU[["647_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.647, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.647, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_647_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.647, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_647_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.647, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_647_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["647_LU"]])


#mag 648
pairwise_p_values_648_LU <- c('0.1041','0.0918','0.2206')


# Create boxplot
plot_list_LU[["648_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.648, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.648, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_648_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.648, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_648_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.648, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_648_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["648_LU"]])

#mag 665
pairwise_p_values_665_LU <- c('0.3932','0.9164','0.3822')


# Create boxplot
plot_list_LU[["665_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.665, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.665, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_665_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.665, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_665_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.665, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_665_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["665_LU"]])

#mag 715
pairwise_p_values_715_LU <- c('0.8852','0.0819','0.0819')


# Create boxplot
plot_list_LU[["715_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.715, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.715, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_715_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.715, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_715_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.715, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_715_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["715_LU"]])


#mag 755
pairwise_p_values_755_LU <- c('0.0359','0.5818','0.0359')


# Create boxplot
plot_list_LU[["755_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.755, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.755, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_755_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.755, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_755_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.755, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_755_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["755_LU"]])


#mag 775
pairwise_p_values_775_LU <- c('0.7087','0.0901','0.0553')


# Create boxplot
plot_list_LU[["775_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.775, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.775, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_775_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.775, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_775_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.775, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_775_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["775_LU"]])

#mag 797
pairwise_p_values_797_LU <- c('0.8499','0.4400','0.4400')


# Create boxplot
plot_list_LU[["797_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.797, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.797, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_797_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.797, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_797_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.797, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_797_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["797_LU"]])

#mag 818
pairwise_p_values_818_LU <- c('0.6525','0.1012','0.0541')


# Create boxplot
plot_list_LU[["818_LU"]] <- ggplot(combined_df, aes(x = Land_use, y = metabat_bins.818, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.818, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_818_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.818, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_818_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.818, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_818_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

print(plot_list_LU[["818_LU"]])

# Arrange the figures in a 2 by 5 grid
grid_row_height <- 20  # Adjust as needed

# Arrange all figures in a 5 by 2 grid with specified height for each row
grid_plot_LU <- do.call(grid.arrange, c(grobs = plot_list_LU, ncol = 2, nrow = 10))
grid_plot_LU$heights <- unit(rep(grid_row_height, 10), "cm")

# Save the grid as an image (e.g., PNG) with high resolution
ggsave("MAG_LME_grid_image_LU.png", plot = grid_plot_LU, width = 30, height = 100, units = "in", dpi = 300, limitsize = FALSE)


#############boxplot for MAG and season

combined_df$Season <- factor(combined_df$Season, levels = c("Autumn", "Spring", "Summer"))
season_colors <- c("Autumn" = "black", "Spring" = "red", "Summer" = "green")

plot_list_S <- list()


#MAG 14
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_14_S <- c("0.5148", "0.3852", "0.5148")


# Create boxplot
plot_list_S[["14_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.14, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.14, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_14_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.14, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_14_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.14, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_14_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_S[["14_S"]])

#MAG 17

pairwise_p_values_17_S <- c('0.9181','0.9181','0.9181')


# Create boxplot
plot_list_S[["17_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.17, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.17, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_17_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.17, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_17_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.17, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_17_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_S[["17_S"]])


#mag 42

pairwise_p_values_42_S <- c('<0.0001','<0.0001','0.1064')


# Create boxplot
plot_list_S[["42_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.42, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.42, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_42_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.42, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_42_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.42, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_42_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_S[["42_S"]])

#mag 47
pairwise_p_values_47_S <- c('<0.0001','<0.0001','0.0527')


# Create boxplot
plot_list_S[["47_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.47, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.47, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_47_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.47, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_47_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.47, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_47_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_S[["47_S"]])


#mag 166

pairwise_p_values_166_S <- c('0.2532','0.2586','0.0320')


# Create boxplot
plot_list_S[["166_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.166, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.166, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_166_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.166, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_166_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.166, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_166_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_LU[["166_S"]])

#mag 209
pairwise_p_values_209_S <- c('0.0452','0.0609','0.7611')


# Create boxplot
plot_list_S[["209_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.209, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.209, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_209_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.209, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_209_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.209, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_209_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_LU[["209_S"]])


#mag 240
pairwise_p_values_240_S <- c('0.7414','0.7414','0.8746')


# Create boxplot
plot_list_S[["240_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.240, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.240, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_240_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.240, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_240_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.240, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_240_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_S[["240_S"]])

#mag 283
pairwise_p_values_283_S <- c('0.7822','0.7822','0.7822')


# Create boxplot
plot_list_S[["283_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.283, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.283, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_283_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.283, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_283_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.283, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_283_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_S[["283_S"]])

#mag 369
pairwise_p_values_369_S <- c('0.0761','0.4379','0.2098')


# Create boxplot
plot_list_S[["369_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.369, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.369, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_369_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.369, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_369_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.369, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_369_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_S[["369_S"]])


#mag 530
pairwise_p_values_530_S <- c('0.0249','<0.0001','0.0025')


# Create boxplot
plot_list_S[["530_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.530, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.530, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_530_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.530, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_530_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.530, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_530_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_S[["530_S"]])


#mag 616
pairwise_p_values_616_S <- c('0.5815','0.0260','0.0124')


# Create boxplot
plot_list_S[["616_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.616, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.616, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_616_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.616, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_616_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.616, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_616_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_S[["616_S"]])


#mag 647

pairwise_p_values_647_S <- c('0.9839','0.9839','0.9839')


# Create boxplot
plot_list_S[["647_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.647, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.647, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_647_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.647, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_647_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.647, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_647_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_S[["647_S"]])


#mag 648
pairwise_p_values_648_S <- c('0.5110','0.5110','0.5110')


# Create boxplot
plot_list_S[["648_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.648, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.648, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_648_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.648, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_648_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.648, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_648_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_S[["648_S"]])

#mag 665
pairwise_p_values_665_S <- c('0.0116','<0.0001','0.0116')


# Create boxplot
plot_list_S[["665_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.665, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.665, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_665_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.665, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_665_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.665, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_665_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_S[["665_S"]])

#mag 715
pairwise_p_values_715_S <- c('0.0615','0.0001','<0.0001')


# Create boxplot
plot_list_S[["715_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.715, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.715, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_715_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.715, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_715_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.715, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_715_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_S[["715_S"]])


#mag 755
pairwise_p_values_755_S <- c('0.0003','<0.0001','0.1654')


# Create boxplot
plot_list_S[["755_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.755, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.755, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_755_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.755, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_755_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.755, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_755_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_S[["755_S"]])


#mag 775
pairwise_p_values_775_S <- c('0.6818','<0.0001','<0.0001')


# Create boxplot
plot_list_S[["775_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.775, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.775, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_775_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.775, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_775_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.775, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_775_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_S[["775_S"]])

#mag 797
pairwise_p_values_797_S <- c('0.7454','0.7454','0.7454')


# Create boxplot
plot_list_S[["797_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.797, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.797, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_797_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.797, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_797_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.797, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_797_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_S[["797_S"]])

#mag 818
pairwise_p_values_818_S <- c('0.1078','<0.0001','0.0001')


# Create boxplot
plot_list_S[["818_S"]] <- ggplot(combined_df, aes(x = Season, y = metabat_bins.818, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.818, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_818_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.818, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_818_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.818, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_818_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors) 

print(plot_list_S[["818_S"]])

# Arrange the figures in a 2 by 5 grid
grid_row_height <- 20  # Adjust as needed

# Arrange all figures in a 5 by 2 grid with specified height for each row
grid_plot_S <- do.call(grid.arrange, c(grobs = plot_list_S, ncol = 2, nrow = 10))
grid_plot_S$heights <- unit(rep(grid_row_height, 10), "cm")

# Save the grid as an image (e.g., PNG) with high resolution
ggsave("MAG_LME_grid_image_S.png", plot = grid_plot_S, width = 30, height = 100, units = "in", dpi = 300, limitsize = FALSE)

##################making MAG vs Interaction effects box plots
#put in grid 
plot_list_IE <- list()
# Ensure Land_use and Season are factors
combined_df$Land_use <- factor(combined_df$Land_use, levels = c("Agri_ditch", "Mixed", "Forested"))
combined_df$Season <- factor(combined_df$Season, levels = c("Autumn", "Spring", "Summer"))

#pre-calculated pairwise pvalue from LME
pairwise_p_values_14_IE <- c("1",
                             "0.4405",
                             "1",
                             "0.2822",
                             "0.4710",
                             "1",
                             "1",
                             "1",
                             "1")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_17_IE <- c("1",
                             "1",
                             "1",
                             "1",
                             "1",
                             "1",
                             "0.8707",
                             "1",
                             "1")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_42_IE <- c("0.0479",
                                  "0.0186",
                                  "0.0052",
                                  "0.0002",
                                  "0.0107",
                                  "0.0021",
                                  "0.0001",
                                  "0.0119",
                                  "0.0031")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_47_IE <- c("1",
                                  "0.1908",
                                  "0.1890",
                                  "0.0041",
                                  "0.1530",
                                  "0.0280",
                                  "0.0010",
                                  "0.2750",
                                  "0.0778")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_166_IE <- c("1",
                              "1",
                              "1",
                              "0.8731",
                              "1",
                              "0.2622",
                              "0.6132",
                              "1",
                              "1")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_209_IE <- c("1",
                                  "0.8922",
                                  "1",
                                  "1",
                                  "1",
                                  "1",
                                  "1",
                                  "1",
                                  "1")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_240_IE <- c("0.1799",
                                "1",
                                "1",
                                "1",
                                "1",
                                "1",
                                "1",
                                "1",
                                "1")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_283_IE <- c("0.3246",
                                    "1",
                                    "0.2616",
                                    "0.4207",
                                    "1",
                                    "0.5184",
                                    "0.1510",
                                    "1",
                                    "0.7823")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_369_IE <- c("1",
                                 "0.7213",
                                 "1",
                                 "1",
                                 "1",
                                 "0.5120",
                                 "0.7738",
                                 "0.3850",
                                 "0.2125")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_530_IE <- c("0.0127",
                                     "0.0835",
                                     "0.0073",
                                     "0.0005",
                                     "0.0326",
                                     "0.0033",
                                     "0.0007",
                                     "0.0781",
                                     "0.0105")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_616_IE <- c("1",
                              "0.0321",
                              "1",
                              "1",
                              "1",
                              "1",
                              "1",
                              "0.4812",
                              "1")

# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_647_IE <- c("1",
                              "0.8718",
                              "1",
                              "1",
                              "0.1247",
                              "1",
                              "1",
                              "0.3863",
                              "1")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_648_IE <- c("1",
                              "0.2611",
                              "1",
                              "0.3341",
                              "1",
                              "1",
                              "0.7766",
                              "0.3503",
                              "1")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_665_IE <- c("0.0005",
                              "0.0235",
                              "0.0037",
                              "6.16e-5",
                              "0.0156",
                              "0.0023",
                              "9.18e-5",
                              "0.0139",
                              "0.0051")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_715_IE <- c("0.0026",
                              "0.1011",
                              "0.1781",
                              "0.0311",
                              "0.4522",
                              "1",
                              "0.0047",
                              "0.1316",
                              "0.0692")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_755_IE <- c("1",
                              "0.0162",
                              "0.0525",
                              "0.0013",
                              "0.0451",
                              "0.0083",
                              "0.0010",
                              "0.0410",
                              "0.0133")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_775_IE <- c("0.0007",
                              "0.1738",
                              "0.2657",
                              "0.3490",
                              "0.7783",
                              "1",
                              "0.02133",
                              "0.1936",
                              "0.1434")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_797_IE <- c("1",
                              "1",
                              "1",
                              "1",
                              "1",
                              "0.8034",
                              "1",
                              "1",
                              "0.9626")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_818_IE <- c("0.0003",
                              "0.0929",
                              "0.0879",
                              "0.0152",
                              "0.2815",
                              "0.5283",
                              "0.0011",
                              "0.0665",
                              "0.0355")

#making the boxplots for interaction effects
#MAG 14
plot_list_IE[["14_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.14, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.14, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_14_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.14, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_14_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.14, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_14_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.14, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_14_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.14, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_14_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.14, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_14_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.14, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_14_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.14, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_14_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.14, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_14_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 14") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#MAG 17
plot_list_IE[["17_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.17, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.17, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_17_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.17, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_17_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.17, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_17_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.17, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_17_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.17, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_17_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.17, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_17_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.17, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_17_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.17, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_17_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.17, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_17_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 17") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#MAG 42
plot_list_IE[["42_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.42, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.42, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_42_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.42, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_42_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.42, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_42_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.42, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_42_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.42, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_42_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.42, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_42_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.42, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_42_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.42, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_42_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.42, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_42_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 42") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#MAG 47
plot_list_IE[["47_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.47, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.47, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_47_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.47, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_47_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.47, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_47_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.47, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_47_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.47, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_47_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.47, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_47_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.47, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_47_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.47, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_47_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.47, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_47_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 47") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#MAG 166
plot_list_IE[["166_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.166, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.166, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_166_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.166, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_166_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.166, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_166_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.166, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_166_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.166, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_166_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.166, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_166_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.166, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_166_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.166, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_166_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.166, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_166_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 166") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#MAG 209
plot_list_IE[["209_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.209, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.209, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_209_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.209, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_209_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.209, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_209_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.209, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_209_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.209, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_209_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.209, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_209_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.209, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_209_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.209, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_209_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.209, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_209_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 209") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))


#MAG 240
plot_list_IE[["240_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.240, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.240, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_240_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.240, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_240_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.240, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_240_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.240, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_240_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.240, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_240_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.240, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_240_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.240, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_240_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.240, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_240_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.240, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_240_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 240") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#MAG 283
plot_list_IE[["283_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.283, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.283, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_283_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.283, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_283_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.283, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_283_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.283, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_283_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.283, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_283_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.283, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_283_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.283, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_283_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.283, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_283_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.283, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_283_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 283") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#MAG 369
plot_list_IE[["369_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.369, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.369, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_369_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.369, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_369_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.369, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_369_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.369, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_369_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.369, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_369_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.369, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_369_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.369, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_369_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.369, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_369_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.369, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_369_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 369") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#MAG 530
plot_list_IE[["530_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.530, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.530, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_530_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.530, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_530_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.530, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_530_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.530, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_530_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.530, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_530_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.530, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_530_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.530, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_530_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.530, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_530_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.530, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_530_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 530") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#MAG 616
plot_list_IE[["616_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.616, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.616, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_616_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.616, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_616_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.616, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_616_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.616, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_616_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.616, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_616_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.616, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_616_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.616, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_616_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.616, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_616_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.616, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_616_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 616") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#MAG 647
plot_list_IE[["647_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.647, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.647, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_647_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.647, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_647_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.647, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_647_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.647, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_647_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.647, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_647_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.647, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_647_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.647, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_647_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.647, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_647_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.647, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_647_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 647") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#MAG 648
plot_list_IE[["648_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.648, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.648, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_648_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.648, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_648_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.648, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_648_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.648, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_648_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.648, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_648_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.648, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_648_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.648, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_648_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.648, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_648_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.648, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_648_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 648") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#MAG 665
plot_list_IE[["665_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.665, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.665, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_665_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.665, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_665_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.665, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_665_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.665, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_665_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.665, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_665_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.665, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_665_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.665, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_665_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.665, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_665_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.665, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_665_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 665") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#MAG 715
plot_list_IE[["715_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.715, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.715, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_715_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.715, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_715_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.715, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_715_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.715, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_715_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.715, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_715_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.715, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_715_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.715, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_715_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.715, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_715_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.715, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_715_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 715") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#MAG 755
plot_list_IE[["755_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.755, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.755, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_755_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.755, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_755_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.755, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_755_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.755, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_755_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.755, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_755_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.755, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_755_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.755, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_755_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.755, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_755_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.755, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_755_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 755") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#MAG 775
plot_list_IE[["775_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.775, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.775, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_775_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.775, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_775_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.775, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_775_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.775, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_775_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.775, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_775_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.775, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_775_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.775, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_775_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.775, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_775_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.775, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_775_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 775") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#MAG 797
plot_list_IE[["797_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.797, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.797, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_797_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.797, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_797_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.797, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_797_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.797, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_797_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.797, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_797_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.797, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_797_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.797, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_797_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.797, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_797_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.797, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_797_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 797") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#MAG 818
plot_list_IE[["818_IE"]] <- ggplot(combined_df, aes(x = interaction(Land_use, Season), y = metabat_bins.818, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(combined_df$metabat_bins.818, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_818_IE[1])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(combined_df$metabat_bins.818, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_818_IE[2])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(combined_df$metabat_bins.818, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_818_IE[3])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(combined_df$metabat_bins.818, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_818_IE[4])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(combined_df$metabat_bins.818, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_818_IE[5])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(combined_df$metabat_bins.818, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_818_IE[6])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(combined_df$metabat_bins.818, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_818_IE[7])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(combined_df$metabat_bins.818, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_818_IE[8])),
            color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(combined_df$metabat_bins.818, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_818_IE[9])),
            color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "MAG 818") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

# Arrange the figures in a 2 by 5 grid
grid_row_height <- 20  # Adjust as needed

# Arrange all figures in a 5 by 2 grid with specified height for each row
grid_plot_IE <- do.call(grid.arrange, c(grobs = plot_list_IE, ncol = 2, nrow = 10))
grid_plot_IE$heights <- unit(rep(grid_row_height, 10), "cm")

# Save the grid as an image (e.g., PNG) with high resolution
ggsave("MAG_LME_grid_image_IE.png", plot = grid_plot_IE, width = 30, height = 100, units = "in", dpi = 300, limitsize = FALSE)

