# Getting the normalized and transformed relative abundances for each taxonomic level, doing Linear Mixed Models and PLS-DA for samples and metadata
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1

#set your working directory were everything will be stored and saved
setwd("C:/Users/Sandersonh/OneDrive - AGR-AGR/Desktop/Wen Chen/Project_1_landuse_AMR_virulence/Statistics_Analysis")


#############################Importing datasets for taxonomy########################################################

# Read data from your TSV files from snakemake pipeline
genus_data_nn <- read.delim("genus.tsv", header = TRUE, stringsAsFactors = FALSE)
#species_data_nn <- read.delim("species.tsv", header = TRUE, stringsAsFactors = FALSE)
family_data_nn <- read.delim("family.tsv", header = TRUE, stringsAsFactors = FALSE)
order_data_nn <- read.delim("order.tsv", header = TRUE, stringsAsFactors = FALSE)

library(tidyr)
library(dplyr)
#calculate relative abundances so the rows add to 1 in crosstable
genus_data_relabund_nn <- genus_data_nn %>%
  group_by(SAMPLE_ID) %>%
  mutate(RelativeAbundance = READS_COUNT / sum(READS_COUNT))
family_data_relabund_nn <- family_data_nn %>%
  group_by(SAMPLE_ID) %>%
  mutate(RelativeAbundance = READS_COUNT / sum(READS_COUNT))
order_data_relabund_nn <- order_data_nn %>%
  group_by(SAMPLE_ID) %>%
  mutate(RelativeAbundance = READS_COUNT / sum(READS_COUNT))

#create crosstables so rows will add to 1

cross_table_genus_nn <- genus_data_relabund_nn %>%
  select(SAMPLE_ID, TAXONOMY, RelativeAbundance) %>%
  spread(TAXONOMY, RelativeAbundance, fill = 0)

cross_table_family_nn <- family_data_relabund_nn %>%
  select(SAMPLE_ID, TAXONOMY, RelativeAbundance) %>%
  spread(TAXONOMY, RelativeAbundance, fill = 0)

cross_table_order_nn <- order_data_relabund_nn %>%
  select(SAMPLE_ID, TAXONOMY, RelativeAbundance) %>%
  spread(TAXONOMY, RelativeAbundance, fill = 0)


#CLR transformations of taxonomic data (Genus, Family, Order)

# Install and load the compositions package if not already installed
#if (!requireNamespace("compositions", quietly = TRUE)) {
#  install.packages("compositions")
#}
library(compositions)
sample_names <- cross_table_genus_nn$SAMPLE_ID
taxa_data_genus <- cross_table_genus_nn[, -which(colnames(cross_table_genus_nn) == "SAMPLE_ID")]
#CLR tranformation for Genus taxonomic data to generate a dataframe
# Perform CLR transformation
clr_data_genus <- clr(taxa_data_genus)

# Convert CLR-transformed data to a data frame
clr_data_genus_df <- as.data.frame(clr_data_genus)

# Add sample names as a column
clr_data_genus_df$SAMPLE_ID <- sample_names

# Reorder columns (optional, adjust as needed)
clr_data_genus_df <- clr_data_genus_df[, c("SAMPLE_ID", colnames(clr_data_genus_df)[1:(ncol(clr_data_genus_df)-1)])]

# Display the CLR-transformed data as a data frame
print(clr_data_genus_df)

#CLR transforamtion for Family taxonomic data to generate a dataframe
sample_names <- cross_table_family_nn$SAMPLE_ID
taxa_data_family <- cross_table_family_nn[, -which(colnames(cross_table_family_nn) == "SAMPLE_ID")]

# Perform CLR transformation
clr_data_family <- clr(taxa_data_family)

# Convert CLR-transformed data to a data frame
clr_data_family_df <- as.data.frame(clr_data_family)

# Add sample names as a column
clr_data_family_df$SAMPLE_ID <- sample_names

# Reorder columns (optional, adjust as needed)
clr_data_family_df <- clr_data_family_df[, c("SAMPLE_ID", colnames(clr_data_family_df)[1:(ncol(clr_data_family_df)-1)])]

# Display the CLR-transformed data as a data frame
print(clr_data_family_df)

#CLR transformation for order taxonomic data to generate a dataframe
sample_names <- cross_table_order_nn$SAMPLE_ID
taxa_data_order <- cross_table_order_nn[, -which(colnames(cross_table_order_nn) == "SAMPLE_ID")]

# Perform CLR transformation
clr_data_order <- clr(taxa_data_order)

# Convert CLR-transformed data to a data frame
clr_data_order_df <- as.data.frame(clr_data_order)

# Add sample names as a column
clr_data_order_df$SAMPLE_ID <- sample_names

# Reorder columns (optional, adjust as needed)
clr_data_order_df <- clr_data_order_df[, c("SAMPLE_ID", colnames(clr_data_order_df)[1:(ncol(clr_data_order_df)-1)])]

# Display the CLR-transformed data as a data frame
print(clr_data_order_df)


############################################Generating large table with all data for statistical analysis#################################################
#install and load packages
#if (!requireNamespace("nlme", quietly = TRUE)) {
#  install.packages("nlme")
#}
#if (!requireNamespace("grafify", quietly = TRUE)) {
#  install.packages("grafify")
#}
#install.packages(c("nlme", "mgcv", "ape"), dependencies = TRUE)
#detach("mgcv", unload = TRUE)
#detach("ape", unload = TRUE)
library(mgcv)
library(ape)
library(nlme)
library(grafify)


#creating a table with all the data we will use, taxonomic data
merged_table <- merge(clr_data_genus_df, clr_data_family_df , by = "SAMPLE_ID", all = TRUE)
merged_table <- merge(merged_table,clr_data_order_df, by = "SAMPLE_ID", all = TRUE)


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


AMR_VF <-read.csv("filtered_AMRVF_crosstab.csv", header = TRUE, stringsAsFactors = FALSE)
metadata <- read.csv("Metadata_clean_AMRVFMG_project.csv", header = TRUE, stringsAsFactors = FALSE) 

metadata <- metadata %>%
  rename(SAMPLE_ID = Sample_ID_Genome)


#merge taxonomy, AMRVF and metadata tables together
merged_table_all_clr <- merge(metadata, AMR_VF, by = "SAMPLE_ID", all = TRUE)
merged_table_all_clr <- merge(merged_table_all_clr, filtered_table, by = "SAMPLE_ID", all = TRUE)

filtered_table_all_clr <- na.omit(merged_table_all_clr)

#Create blocks with SN18 ans SN19 as one block, SN20 and SN21 as one block and the rest of the SAmple_Sites as seperate blocks

filtered_table_all_clr$Block <- ifelse(filtered_table_all_clr$Sample_Site %in% c("SN018", "SN019"), "Block1819",
                          ifelse(filtered_table_all_clr$Sample_Site %in% c("SN020", "SN021"), "Block2021", paste0("Block_", filtered_table_all_clr$Sample_Site)))

colnames(filtered_table_all_clr) <- sub(" ", "_", colnames(filtered_table_all_clr))
colnames(filtered_table_all_clr)[colnames(filtered_table_all_clr) == "Pusillimonas_(ex Stolz et al. 2005)"] <- "Pusillimonas"
colnames(filtered_table_all_clr)[colnames(filtered_table_all_clr) == "Clostridiales_Family XVI. Incertae Sedis"] <- "Clostridiales_Family_XVI_Incertae_Sedis"
colnames(filtered_table_all_clr)[colnames(filtered_table_all_clr) == "Clostridiales_Family XVII. Incertae Sedis"] <- "Clostridiales_Family_XVII_Incertae_Sedis"
colnames(filtered_table_all_clr)[colnames(filtered_table_all_clr) == "Eubacteriales_Family XIII. Incertae Sedis"] <- "Eubacteriales_Family_XIII_Incertae_Sedis"
colnames(filtered_table_all_clr)[colnames(filtered_table_all_clr) == "Thermoanaerobacterales_Family III. Incertae Sedis"] <- "Thermoanaerobacterales_Family_III_Incertae_Sedis"
colnames(filtered_table_all_clr)[colnames(filtered_table_all_clr) == "Thermoanaerobacterales_Family IV. Incertae Sedis"] <- "Thermoanaerobacterales_Family_IV_Incertae_Sedis"
colnames(filtered_table_all_clr)[colnames(filtered_table_all_clr) == "Bacteroidetes_Order II. Incertae sedis"] <- "Bacteroidetes_Order_II_Incertae_Sedis"


#write.csv(filtered_table_all_clr, file = "table_lme_metadata_AMRVF_taxa.csv", row.names = FALSE)





#############################################LME##############################################################
#Linear mixed effects models to assess effects of land use types on water physicochemical properties, abundance of taxa at different ranks, presence and absence of AMR/VF
#significant level of p=<0.05
#land use as fixed effect and block+samplingdata as random effects; sampling date as a repeated measure
#pairwise comparison using posthoc_pairwise function in grafify with "tukey adjustment
#partial least squares discriminant analysis (PLS-DA) from DiscriMiner to discriminate SNR sampling sites based on independent metadata data
#discriminer is not available for my R version will need to find another package

#land use as fixed effect and block+samplingdata as random effects; sampling date as a repeated measure
# Exclude variables from response variables
library(nlme)
#install.packages("emmeans")
library(emmeans)
# Assuming your "date" variable is in Date format, convert it to numeric (days since a reference date)
filtered_table_all_clr$date <- as.Date(filtered_table_all_clr$date, format = "%m/%d/%y")
filtered_table_all_clr$date_numeric <- as.numeric(filtered_table_all_clr$date - min(filtered_table_all_clr$date))

AMR_VF_column_names <- colnames(AMR_VF)
exclude_variables <- c("Sample_ID", "Block", "date", "desc_", "Land_use", "sample_type", "Season", "sampleID", "site_yr", "yeardayofyear", "Sample_Site", "Site", "Year", "Week","site", "strahler", "lon", "lat", AMR_VF_column_names)
response_variables <- setdiff(names(filtered_table_all_clr), exclude_variables)

# Open a file connection for writing
file_path <- "lme_contdate_landuse_results.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ Land_use", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = filtered_table_all_clr)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      print(summary(model))
     
       # ANOVA
      cat("ANOVA:\n")
      print(anova(model))
      
      # Pairwise comparisons
      cat("Pairwise Comparison Results:\n")
      pwc <- emmeans(model, specs = "Land_use", data = filtered_table_all_clr)
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


#Season as the fixed variable
file_path <- "lme_contdate_season_results.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ Season", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = filtered_table_all_clr)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      print(summary(model))
      
      # ANOVA
      cat("ANOVA:\n")
      print(anova(model))
      
      # Pairwise comparisons
      cat("Pairwise Comparison Results:\n")
      pwc <- emmeans(model, specs = "Season", data = filtered_table_all_clr)
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

#Week as fixed variable

file_path <- "lme_contdate_week_results.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ Week", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = filtered_table_all_clr)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      print(summary(model))
      
      # ANOVA
      cat("ANOVA:\n")
      print(anova(model))
      
      # Pairwise comparisons
      cat("Pairwise Comparison Results:\n")
      pwc <- emmeans(model, specs = "Week", data = filtered_table_all_clr)
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


#Year as fixed variable
filtered_table_all_clr$Year <- as.factor(filtered_table_all_clr$Year)
file_path <- "lme_contdate_year_results.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ Year", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = filtered_table_all_clr)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      print(summary(model))
      
      # ANOVA
      cat("ANOVA:\n")
      print(anova(model))
      
      # Pairwise comparisons
      cat("Pairwise Comparison Results:\n")
      pwc <- emmeans(model, specs = "Year", data = filtered_table_all_clr)
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

#month as fixed variable
filtered_table_all_clr <- read.csv("table_lme_metadata_AMRVF_taxa_month.csv", header = TRUE, stringsAsFactors = FALSE)

AMR_VF_column_names <- colnames(AMR_VF)
exclude_variables <- c("Sample_ID", "Block", "date", "desc_", "Land_use", "sample_type", "Season", "sampleID", "site_yr", "yeardayofyear", "Sample_Site", "Site", "Year", "Week", "Month", "site", "strahler", "lon", "lat", "data_numeric", AMR_VF_column_names)
response_variables <- setdiff(names(filtered_table_all_clr), exclude_variables)
file_path <- "lme_contdate_month_results.txt"
filtered_table_all_clr$date <- as.Date(filtered_table_all_clr$date, format = "%m/%d/%y")
filtered_table_all_clr$date_numeric <- as.numeric(filtered_table_all_clr$date - min(filtered_table_all_clr$date))

capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ Month", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = filtered_table_all_clr)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      print(summary(model))
      
      # ANOVA
      cat("ANOVA:\n")
      print(anova(model))
      
      # Pairwise comparisons
      cat("Pairwise Comparison Results:\n")
      pwc <- emmeans(model, specs = "Month", data = filtered_table_all_clr)
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

filtered_table_all_clr$Year <- as.factor(filtered_table_all_clr$Year)


#AMR/VFs as fixed variable for lme
#change all AMR/VF columns to be seen as factors
filtered_table_all_clr[AMR_VF_column_names] <- lapply(filtered_table_all_clr[AMR_VF_column_names], factor)


#TETM
file_path <- "lme_contdate_tetM_results.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ TETM", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = filtered_table_all_clr)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      print(summary(model))
      
      # ANOVA
      cat("ANOVA:\n")
      print(anova(model))
      
      # Pairwise comparisons
      cat("Pairwise Comparison Results:\n")
      pwc <- emmeans(model, specs = "TETM", data = filtered_table_all_clr)
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

#LSA

file_path <- "lme_contdate_lsa_results.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ LSA", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = filtered_table_all_clr)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      print(summary(model))
      
      # ANOVA
      cat("ANOVA:\n")
      print(anova(model))
      
      # Pairwise comparisons
      cat("Pairwise Comparison Results:\n")
      pwc <- emmeans(model, specs = "LSA", data = filtered_table_all_clr)
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

#LSAE

file_path <- "lme_contdate_lsae_results.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ LSAE", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = filtered_table_all_clr)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      print(summary(model))
      
      # ANOVA
      cat("ANOVA:\n")
      print(anova(model))
      
      # Pairwise comparisons
      cat("Pairwise Comparison Results:\n")
      pwc <- emmeans(model, specs = "LSAE", data = filtered_table_all_clr)
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

#SODB

file_path <- "lme_contdate_sodb_results.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ SODB", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = filtered_table_all_clr)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      print(summary(model))
      
      # ANOVA
      cat("ANOVA:\n")
      print(anova(model))
      
      # Pairwise comparisons
      cat("Pairwise Comparison Results:\n")
      pwc <- emmeans(model, specs = "SODB", data = filtered_table_all_clr)
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

#VFG010906
file_path <- "lme_contdate_VFG010906_results.txt"
capture.output({
  sink(file_path, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ VFG010906", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = filtered_table_all_clr)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      print(summary(model))
      
      # ANOVA
      cat("ANOVA:\n")
      print(anova(model))
      
      # Pairwise comparisons
      cat("Pairwise Comparison Results:\n")
      pwc <- emmeans(model, specs = "VFG010906", data = filtered_table_all_clr)
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
##############GEtting pvalues for ANOVA and Pairwise Comparisons
# Load necessary libraries
# Load necessary libraries


library(nlme)
library(emmeans)
library(broom)

# Specify the full path to the file
file_path <- "lme_contdate_landuse_results.txt"
file_path2 <-"lme_contdate_landuse_results_pvalues.txt"

# Function to extract response variables from the file
extract_response_variables <- function(file_path) {
  response_variables <- character()
  lines <- readLines(file_path)
  for (line in lines) {
    if (grepl("Response Variable:", line)) {
      response_variable <- gsub("^Response Variable:\\s+", "", line)
      response_variables <- c(response_variables, response_variable)
    }
  }
  return(response_variables)
}

# Get response variables
response_variables <- extract_response_variables(file_path)

# Create empty data frames to store results
anova_df <- data.frame(Response_Variable = character(),
                       ANOVA_PValue = numeric(),
                       stringsAsFactors = FALSE)

pairwise_df <- data.frame(Response_Variable = character(),
                          Comparison = character(),
                          Estimate = numeric(),
                          SE = numeric(),
                          df = numeric(),
                          conf.low = numeric(),
                          conf.high = numeric(),
                          statistic = numeric(),
                          PValue = numeric(),
                          stringsAsFactors = FALSE)

# Capture output to the file
capture.output({
  # Open the file connection
  sink(file_path2, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ Land_use", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = filtered_table_all_clr)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      #print(summary(model))
      
      # ANOVA
      cat("ANOVA:\n")
      anova_result <- anova(model)
      df_anova_result <- data.frame(anova_result)
      df_anova_result["Land_use","p.value"]
    
      
      # Extract ANOVA p-value and save it to the file and data frame
      anova_p_value <- anova_result$"p-value"[2]
      cat("ANOVA p-value: ", anova_p_value, "\n")
      anova_df <- rbind(anova_df, data.frame(Response_Variable = response_variable,
                                             ANOVA_PValue = anova_p_value,
                                             stringsAsFactors = FALSE))
      
      # Pairwise comparisons
      cat("Pairwise Comparison Results:\n")
      pwc <- emmeans(model, specs = "Land_use", data = filtered_table_all_clr)
      pairwise_results <- tidy(pwc, conf.int = TRUE)
      print(pairwise_results)
      
      # Save pairwise results to data frame
      pairwise_df <- rbind(pairwise_df, data.frame(Response_Variable = response_variable,
                                                   Comparison = pairwise_results[,1],
                                                   Estimate = pairwise_results$estimate,
                                                   SE = pairwise_results$std.error,
                                                   df = pairwise_results$df,
                                                   conf.low = pairwise_results$conf.low,
                                                   cong.high = pairwise_results$conf.high,
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
}, file = file_path2, append = TRUE)

# Check the file contents
cat(readLines(file_path), sep = "\n")

# Print data frames
print(anova_df)
print(pairwise_df)

write.csv(anova_df, file = "lme_ANOVA_pvalues_landuse.csv", row.names = FALSE)
write.csv(pairwise_df, file = "lme_PWC_pvalues_landuse.csv", row.names = FALSE)

####Combine the two tables
anova <- read.csv("lme_ANOVA_pvalues_landuse.csv", header = TRUE, stringsAsFactors = FALSE) 
pairwise <- read.csv("lme_PWC_pvalues_landuse.csv", header = TRUE, stringsAsFactors = FALSE) 
merged_table_pvalues <- merge(anova, pairwise, by.x = 'Response_Variable', by.y = 'Response_Variable')
write.csv(merged_table_pvalues, file = "lme_ALL_pvalues_landuse.csv", row.names = FALSE)

library(dplyr)

# Assuming you have your merged_table data frame
filtered_table <- merged_table_pvalues %>%
  filter(ANOVA_PValue < 0.05)

# Display the result
print(filtered_table)
write.csv(filtered_table, file = "lme_ALL_pvalues_landuse_only_significant_ANOVA.csv", row.names = FALSE)


###Season####

# Specify the full path to the file
file_path <- "lme_contdate_Season_results.txt"
file_path2 <-"lme_contdate_Season_results_pvalues.txt"

# Function to extract response variables from the file
extract_response_variables <- function(file_path) {
  response_variables <- character()
  lines <- readLines(file_path)
  for (line in lines) {
    if (grepl("Response Variable:", line)) {
      response_variable <- gsub("^Response Variable:\\s+", "", line)
      response_variables <- c(response_variables, response_variable)
    }
  }
  return(response_variables)
}

# Get response variables
response_variables <- extract_response_variables(file_path)

# Create empty data frames to store results
anova_df <- data.frame(Response_Variable = character(),
                       ANOVA_PValue = numeric(),
                       stringsAsFactors = FALSE)

pairwise_df <- data.frame(Response_Variable = character(),
                          Comparison = character(),
                          Estimate = numeric(),
                          SE = numeric(),
                          df = numeric(),
                          conf.low = numeric(),
                          conf.high = numeric(),
                          statistic = numeric(),
                          PValue = numeric(),
                          stringsAsFactors = FALSE)

# Capture output to the file
capture.output({
  # Open the file connection
  sink(file_path2, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ Season", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = filtered_table_all_clr)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      
      # ANOVA
      cat("ANOVA:\n")
      anova_result <- anova(model)
      df_anova_result <- data.frame(anova_result)
      df_anova_result["Season","p.value"]
      
      
      # Extract ANOVA p-value and save it to the file and data frame
      anova_p_value <- anova_result$"p-value"[2]
      cat("ANOVA p-value: ", anova_p_value, "\n")
      anova_df <- rbind(anova_df, data.frame(Response_Variable = response_variable,
                                             ANOVA_PValue = anova_p_value,
                                             stringsAsFactors = FALSE))
      
      # Pairwise comparisons
      cat("Pairwise Comparison Results:\n")
      pwc <- emmeans(model, specs = "Season", data = filtered_table_all_clr)
      pairwise_results <- tidy(pwc, conf.int = TRUE)
      print(pairwise_results)
      
      # Save pairwise results to data frame
      pairwise_df <- rbind(pairwise_df, data.frame(Response_Variable = response_variable,
                                                   Comparison = pairwise_results[, 1],
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
}, file = file_path2, append = TRUE)

# Check the file contents
cat(readLines(file_path2), sep = "\n")

# Print data frames
print(anova_df)
print(pairwise_df)

write.csv(anova_df, file = "lme_ANOVA_pvalues_season.csv", row.names = FALSE)
write.csv(pairwise_df, file = "lme_PWC_pvalues_season.csv", row.names = FALSE)
####Combine the two tables
anova_S <- read.csv("lme_ANOVA_pvalues_season.csv", header = TRUE, stringsAsFactors = FALSE) 
pairwise_S <- read.csv("lme_PWC_pvalues_season.csv", header = TRUE, stringsAsFactors = FALSE) 
merged_table_pvalues_S <- merge(anova_S, pairwise_S, by.x = 'Response_Variable', by.y = 'Response_Variable')
write.csv(merged_table_pvalues_S, file = "lme_ALL_pvalues_season.csv", row.names = FALSE)

library(dplyr)

# Assuming you have your merged_table data frame
filtered_table_S <- merged_table_pvalues_S %>%
  filter(ANOVA_PValue < 0.05)

# Display the result
print(filtered_table_S)
write.csv(filtered_table_S, file = "lme_ALL_pvalues_season_only_significant_ANOVA.csv", row.names = FALSE)

##################Interaction Effects--Land_Use and Season###########################

# Specify the full path to the file
file_path <- "lme_contdate_season_results.txt"
file_path2 <-"lme_contdate_seasonlanduse_results_pvalues.txt"

# Function to extract response variables from the file
extract_response_variables <- function(file_path) {
  response_variables <- character()
  lines <- readLines(file_path)
  for (line in lines) {
    if (grepl("Response Variable:", line)) {
      response_variable <- gsub("^Response Variable:\\s+", "", line)
      response_variables <- c(response_variables, response_variable)
    }
  }
  return(response_variables)
}

# Get response variables
response_variables <- extract_response_variables(file_path)

# Create empty data frames to store results
anova_df <- data.frame(Response_Variable = character(),
                       Land_Use_PValue = numeric(),
                       Season_Pvalue = numeric(),
                       Interaction_Pvalue = numeric(),
                       stringsAsFactors = FALSE)

pairwise_df <- data.frame(Response_Variable = character(),
                          Comparison1 = character(),
                          Comparison2 = character(),
                          Estimate = numeric(),
                          SE = numeric(),
                          df = numeric(),
                          conf.low = numeric(),
                          conf.high = numeric(),
                          statistic = numeric(),
                          PValue = numeric(),
                          stringsAsFactors = FALSE)

# Capture output to the file
capture.output({
  # Open the file connection
  sink(file_path2, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with lme, including the interaction term
    formula_str <- paste(response_variable, "~ Land_use * Season", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(fixed = formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = filtered_table_all_clr)
      
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
}, file = file_path2, append = TRUE)

# Check the file contents
cat(readLines(file_path2), sep = "\n")

# Print data frames
print(anova_df)
print(pairwise_df)

write.csv(anova_df, file = "lme_ANOVA_pvalues_seasonlanduse.csv", row.names = FALSE)
write.csv(pairwise_df, file = "lme_PWC_pvalues_seasonlanduse.csv", row.names = FALSE)
####Combine the two tables
anova_S <- read.csv("lme_ANOVA_pvalues_seasonlanduse.csv", header = TRUE, stringsAsFactors = FALSE) 
pairwise_S <- read.csv("lme_PWC_pvalues_seasonlanduse.csv", header = TRUE, stringsAsFactors = FALSE) 
merged_table_pvalues_S <- merge(anova_S, pairwise_S, by.x = 'Response_Variable', by.y = 'Response_Variable')
write.csv(merged_table_pvalues_S, file = "lme_ALL_pvalues_seasonlanduse.csv", row.names = FALSE)

library(dplyr)

# Assuming you have your merged_table data frame
filtered_table_S <- merged_table_pvalues_S %>%
  filter(Interaction_PValue < 0.05)

# Display the result
print(filtered_table_S)
write.csv(filtered_table_S, file = "lme_ALL_pvalues_seasonlanduse_only_significant_ANOVA.csv", row.names = FALSE)

################################AMR/VF LME results###################################
##tetm
# Specify the full path to the file
file_path <- "lme_contdate_tetM_results.txt"
file_path2 <-"lme_contdate_tetM_results_pvalues.txt"

# Function to extract response variables from the file
extract_response_variables <- function(file_path) {
  response_variables <- character()
  lines <- readLines(file_path)
  for (line in lines) {
    if (grepl("Response Variable:", line)) {
      response_variable <- gsub("^Response Variable:\\s+", "", line)
      response_variables <- c(response_variables, response_variable)
    }
  }
  return(response_variables)
}

# Get response variables
response_variables <- extract_response_variables(file_path)

# Create empty data frames to store results
anova_df <- data.frame(Response_Variable = character(),
                       ANOVA_PValue = numeric(),
                       stringsAsFactors = FALSE)

#pairwise_df <- data.frame(Response_Variable = character(),
#                          Comparison = character(),
#                          Estimate = numeric(),
#                          SE = numeric(),
#                          df = numeric(),
#                          conf.low = numeric(),
#                          conf.high = numeric(),
#                          statistic = numeric(),
#                          PValue = numeric(),
#                          stringsAsFactors = FALSE)
#
# Capture output to the file
capture.output({
  # Open the file connection
  sink(file_path2, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ TETM", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = filtered_table_all_clr)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      
      # ANOVA
      cat("ANOVA:\n")
      anova_result <- anova(model)
      df_anova_result <- data.frame(anova_result)
      df_anova_result["TETM","p.value"]
      
      
      # Extract ANOVA p-value and save it to the file and data frame
      anova_p_value <- anova_result$"p-value"[2]
      cat("ANOVA p-value: ", anova_p_value, "\n")
      anova_df <- rbind(anova_df, data.frame(Response_Variable = response_variable,
                                             ANOVA_PValue = anova_p_value,
                                             stringsAsFactors = FALSE))
      
      # Pairwise comparisons
 #     cat("Pairwise Comparison Results:\n")
 #     pwc <- emmeans(model, specs = "TETM", data = filtered_table_all_clr)
 #     pairwise_results <- tidy(pwc, conf.int = TRUE)
 #     print(pairwise_results)
      
      # Save pairwise results to data frame
 #     pairwise_df <- rbind(pairwise_df, data.frame(Response_Variable = response_variable,
 #                                                  Comparison = pairwise_results[, 1],
 #                                                  Estimate = pairwise_results$estimate,
 #                                                  SE = pairwise_results$std.error,
 #                                                  df = pairwise_results$df,
  #                                                 conf.low = pairwise_results$conf.low,
 #                                                  conf.high = pairwise_results$conf.high,
 #                                                  statistic = pairwise_results$statistic,
  #                                                 PValue = p.adjust(pairwise_results$p.value, method = "bonferroni"),
  #                                                 stringsAsFactors = FALSE))
      cat("\n")
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
    # Close the file connection
    sink()
  }, file = file_path2, append = TRUE)
# Check the file contents
cat(readLines(file_path2), sep = "\n")

# Print data frames
print(anova_df)
#print(pairwise_df)


write.csv(anova_df, file = "lme_ANOVA_pvalues_tetM.csv", row.names = FALSE)
#write.csv(pairwise_df, file = "lme_PWC_pvalues_tetM.csv", row.names = FALSE)

####Combine the two tables
anova <- read.csv("lme_ANOVA_pvalues_tetM.csv", header = TRUE, stringsAsFactors = FALSE) 
#pairwise <- read.csv("lme_PWC_pvalues_tetM.csv", header = TRUE, stringsAsFactors = FALSE) 
#merged_table_pvalues <- merge(anova, pairwise, by.x = 'Response_Variable', by.y = 'Response_Variable')
#write.csv(merged_table_pvalues, file = "lme_ALL_pvalues_tetM.csv", row.names = FALSE)

library(dplyr)

# Assuming you have your merged_table data frame
filtered_table <- anova %>%
  filter(ANOVA_PValue < 0.05)

# Display the result
print(filtered_table)
write.csv(filtered_table, file = "lme_ALL_pvalues_tetM_only_significant_ANOVA.csv", row.names = FALSE)
##################LSA###################################
# Specify the full path to the file
file_path <- "lme_contdate_lsa_results.txt"
file_path2 <-"lme_contdate_lsa_results_pvalues.txt"

# Function to extract response variables from the file
extract_response_variables <- function(file_path) {
  response_variables <- character()
  lines <- readLines(file_path)
  for (line in lines) {
    if (grepl("Response Variable:", line)) {
      response_variable <- gsub("^Response Variable:\\s+", "", line)
      response_variables <- c(response_variables, response_variable)
    }
  }
  return(response_variables)
}

# Get response variables
response_variables <- extract_response_variables(file_path)

# Create empty data frames to store results
anova_df <- data.frame(Response_Variable = character(),
                       ANOVA_PValue = numeric(),
                       stringsAsFactors = FALSE)

#pairwise_df <- data.frame(Response_Variable = character(),
#                          Comparison = character(),
#                          Estimate = numeric(),
#                          SE = numeric(),
#                          df = numeric(),
#                          conf.low = numeric(),
#                          conf.high = numeric(),
#                          statistic = numeric(),
#                          PValue = numeric(),
#                          stringsAsFactors = FALSE)

# Capture output to the file
capture.output({
  # Open the file connection
  sink(file_path2, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ LSA", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = filtered_table_all_clr)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      
      # ANOVA
      cat("ANOVA:\n")
      anova_result <- anova(model)
      df_anova_result <- data.frame(anova_result)
      df_anova_result["LSA","p.value"]
      
      
      # Extract ANOVA p-value and save it to the file and data frame
      anova_p_value <- anova_result$"p-value"[2]
      cat("ANOVA p-value: ", anova_p_value, "\n")
      anova_df <- rbind(anova_df, data.frame(Response_Variable = response_variable,
                                             ANOVA_PValue = anova_p_value,
                                             stringsAsFactors = FALSE))
      
      # Pairwise comparisons
 #     cat("Pairwise Comparison Results:\n")
 #     pwc <- emmeans(model, specs = "TETM", data = filtered_table_all_clr)
 #     pairwise_results <- tidy(pwc, conf.int = TRUE)
 #     print(pairwise_results)
      
      # Save pairwise results to data frame
 #     pairwise_df <- rbind(pairwise_df, data.frame(Response_Variable = response_variable,
 #                                                  Comparison = pairwise_results[, 1],
 #                                                  Estimate = pairwise_results$estimate,
 #                                                  SE = pairwise_results$std.error,
 #                                                  df = pairwise_results$df,
 #                                                  conf.low = pairwise_results$conf.low,
 #                                                  conf.high = pairwise_results$conf.high,
 #                                                  statistic = pairwise_results$statistic,
  #                                                 PValue = p.adjust(pairwise_results$p.value, method = "bonferroni"),
  #                                                 stringsAsFactors = FALSE))
      cat("\n")
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
  # Close the file connection
  sink()
}, file = file_path2, append = TRUE)
# Check the file contents
cat(readLines(file_path2), sep = "\n")

# Print data frames
print(anova_df)
#print(pairwise_df)


write.csv(anova_df, file = "lme_ANOVA_pvalues_lsa.csv", row.names = FALSE)
#write.csv(pairwise_df, file = "lme_PWC_pvalues_lsa.csv", row.names = FALSE)

####Combine the two tables
anova <- read.csv("lme_ANOVA_pvalues_lsa.csv", header = TRUE, stringsAsFactors = FALSE) 
#pairwise <- read.csv("lme_PWC_pvalues_lsa.csv", header = TRUE, stringsAsFactors = FALSE) 
#merged_table_pvalues <- merge(anova, pairwise, by.x = 'Response_Variable', by.y = 'Response_Variable')
#write.csv(merged_table_pvalues, file = "lme_ALL_pvalues_lsa.csv", row.names = FALSE)

library(dplyr)

# Assuming you have your merged_table data frame
filtered_table <- anova %>%
  filter(ANOVA_PValue < 0.05)

# Display the result
print(filtered_table)
write.csv(filtered_table, file = "lme_ALL_pvalues_lsa_only_significant_ANOVA.csv", row.names = FALSE)




#####LSAE
# Specify the full path to the file
file_path <- "lme_contdate_lsae_results.txt"
file_path2 <-"lme_contdate_lsae_results_pvalues.txt"

# Function to extract response variables from the file
extract_response_variables <- function(file_path) {
  response_variables <- character()
  lines <- readLines(file_path)
  for (line in lines) {
    if (grepl("Response Variable:", line)) {
      response_variable <- gsub("^Response Variable:\\s+", "", line)
      response_variables <- c(response_variables, response_variable)
    }
  }
  return(response_variables)
}

# Get response variables
response_variables <- extract_response_variables(file_path)

# Create empty data frames to store results
anova_df <- data.frame(Response_Variable = character(),
                       ANOVA_PValue = numeric(),
                       stringsAsFactors = FALSE)

#pairwise_df <- data.frame(Response_Variable = character(),
#                          Comparison = character(),
#                          Estimate = numeric(),
#                          SE = numeric(),
#                          df = numeric(),
#                          conf.low = numeric(),
#                          conf.high = numeric(),
#                          statistic = numeric(),
#                          PValue = numeric(),
#                          stringsAsFactors = FALSE)

# Capture output to the file
capture.output({
  # Open the file connection
  sink(file_path2, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ LSAE", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = filtered_table_all_clr)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      
      # ANOVA
      cat("ANOVA:\n")
      anova_result <- anova(model)
      df_anova_result <- data.frame(anova_result)
      df_anova_result["LSAE","p.value"]
      
      
      # Extract ANOVA p-value and save it to the file and data frame
      anova_p_value <- anova_result$"p-value"[2]
      cat("ANOVA p-value: ", anova_p_value, "\n")
      anova_df <- rbind(anova_df, data.frame(Response_Variable = response_variable,
                                             ANOVA_PValue = anova_p_value,
                                             stringsAsFactors = FALSE))
      
      # Pairwise comparisons
      #     cat("Pairwise Comparison Results:\n")
      #     pwc <- emmeans(model, specs = "TETM", data = filtered_table_all_clr)
      #     pairwise_results <- tidy(pwc, conf.int = TRUE)
      #     print(pairwise_results)
      
      # Save pairwise results to data frame
      #     pairwise_df <- rbind(pairwise_df, data.frame(Response_Variable = response_variable,
      #                                                  Comparison = pairwise_results[, 1],
      #                                                  Estimate = pairwise_results$estimate,
      #                                                  SE = pairwise_results$std.error,
      #                                                  df = pairwise_results$df,
      #                                                  conf.low = pairwise_results$conf.low,
      #                                                  conf.high = pairwise_results$conf.high,
      #                                                  statistic = pairwise_results$statistic,
      #                                                 PValue = p.adjust(pairwise_results$p.value, method = "bonferroni"),
      #                                                 stringsAsFactors = FALSE))
      cat("\n")
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
  # Close the file connection
  sink()
}, file = file_path2, append = TRUE)
# Check the file contents
cat(readLines(file_path2), sep = "\n")

# Print data frames
print(anova_df)
#print(pairwise_df)


write.csv(anova_df, file = "lme_ANOVA_pvalues_lsae.csv", row.names = FALSE)
#write.csv(pairwise_df, file = "lme_PWC_pvalues_lsae.csv", row.names = FALSE)

####Combine the two tables
anova <- read.csv("lme_ANOVA_pvalues_lsae.csv", header = TRUE, stringsAsFactors = FALSE) 
#pairwise <- read.csv("lme_PWC_pvalues_lsae.csv", header = TRUE, stringsAsFactors = FALSE) 
#merged_table_pvalues <- merge(anova, pairwise, by.x = 'Response_Variable', by.y = 'Response_Variable')
#write.csv(merged_table_pvalues, file = "lme_ALL_pvalues_lsae.csv", row.names = FALSE)

library(dplyr)

# Assuming you have your merged_table data frame
filtered_table <- anova %>%
  filter(ANOVA_PValue < 0.05)

# Display the result
print(filtered_table)
write.csv(filtered_table, file = "lme_ALL_pvalues_lsae_only_significant_ANOVA.csv", row.names = FALSE)

######SODB#################
file_path <- "lme_contdate_sodb_results.txt"
file_path2 <-"lme_contdate_sodb_results_pvalues.txt"

# Function to extract response variables from the file
extract_response_variables <- function(file_path) {
  response_variables <- character()
  lines <- readLines(file_path)
  for (line in lines) {
    if (grepl("Response Variable:", line)) {
      response_variable <- gsub("^Response Variable:\\s+", "", line)
      response_variables <- c(response_variables, response_variable)
    }
  }
  return(response_variables)
}

# Get response variables
response_variables <- extract_response_variables(file_path)

# Create empty data frames to store results
anova_df <- data.frame(Response_Variable = character(),
                       ANOVA_PValue = numeric(),
                       stringsAsFactors = FALSE)

#pairwise_df <- data.frame(Response_Variable = character(),
#                          Comparison = character(),
#                          Estimate = numeric(),
#                          SE = numeric(),
#                          df = numeric(),
#                          conf.low = numeric(),
#                          conf.high = numeric(),
#                          statistic = numeric(),
#                          PValue = numeric(),
#                          stringsAsFactors = FALSE)

# Capture output to the file
capture.output({
  # Open the file connection
  sink(file_path2, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ SODB", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = filtered_table_all_clr)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      
      # ANOVA
      cat("ANOVA:\n")
      anova_result <- anova(model)
      df_anova_result <- data.frame(anova_result)
      df_anova_result["SODB","p.value"]
      
      
      # Extract ANOVA p-value and save it to the file and data frame
      anova_p_value <- anova_result$"p-value"[2]
      cat("ANOVA p-value: ", anova_p_value, "\n")
      anova_df <- rbind(anova_df, data.frame(Response_Variable = response_variable,
                                             ANOVA_PValue = anova_p_value,
                                             stringsAsFactors = FALSE))
      
      # Pairwise comparisons
      #     cat("Pairwise Comparison Results:\n")
      #     pwc <- emmeans(model, specs = "TETM", data = filtered_table_all_clr)
      #     pairwise_results <- tidy(pwc, conf.int = TRUE)
      #     print(pairwise_results)
      
      # Save pairwise results to data frame
      #     pairwise_df <- rbind(pairwise_df, data.frame(Response_Variable = response_variable,
      #                                                  Comparison = pairwise_results[, 1],
      #                                                  Estimate = pairwise_results$estimate,
      #                                                  SE = pairwise_results$std.error,
      #                                                  df = pairwise_results$df,
      #                                                  conf.low = pairwise_results$conf.low,
      #                                                  conf.high = pairwise_results$conf.high,
      #                                                  statistic = pairwise_results$statistic,
      #                                                 PValue = p.adjust(pairwise_results$p.value, method = "bonferroni"),
      #                                                 stringsAsFactors = FALSE))
      cat("\n")
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
  # Close the file connection
  sink()
}, file = file_path2, append = TRUE)
# Check the file contents
cat(readLines(file_path2), sep = "\n")

# Print data frames
print(anova_df)
#print(pairwise_df)


write.csv(anova_df, file = "lme_ANOVA_pvalues_sodb.csv", row.names = FALSE)
#write.csv(pairwise_df, file = "lme_PWC_pvalues_sodb.csv", row.names = FALSE)

####Combine the two tables
anova <- read.csv("lme_ANOVA_pvalues_sodb.csv", header = TRUE, stringsAsFactors = FALSE) 
#pairwise <- read.csv("lme_PWC_pvalues_sodb.csv", header = TRUE, stringsAsFactors = FALSE) 
#merged_table_pvalues <- merge(anova, pairwise, by.x = 'Response_Variable', by.y = 'Response_Variable')
#write.csv(merged_table_pvalues, file = "lme_ALL_pvalues_sodb.csv", row.names = FALSE)

library(dplyr)

# Assuming you have your merged_table data frame
filtered_table <- anova %>%
  filter(ANOVA_PValue < 0.05)

# Display the result
print(filtered_table)
write.csv(filtered_table, file = "lme_ALL_pvalues_sodb_only_significant_ANOVA.csv", row.names = FALSE)


###############VFG010906#####################

file_path <- "lme_contdate_VFG010906_results.txt"
file_path2 <-"lme_contdate_VFG010906_results_pvalues.txt"

# Function to extract response variables from the file
extract_response_variables <- function(file_path) {
  response_variables <- character()
  lines <- readLines(file_path)
  for (line in lines) {
    if (grepl("Response Variable:", line)) {
      response_variable <- gsub("^Response Variable:\\s+", "", line)
      response_variables <- c(response_variables, response_variable)
    }
  }
  return(response_variables)
}

# Get response variables
response_variables <- extract_response_variables(file_path)

# Create empty data frames to store results
anova_df <- data.frame(Response_Variable = character(),
                       ANOVA_PValue = numeric(),
                       stringsAsFactors = FALSE)

#pairwise_df <- data.frame(Response_Variable = character(),
#                          Comparison = character(),
#                          Estimate = numeric(),
#                          SE = numeric(),
#                          df = numeric(),
#                          conf.low = numeric(),
#                          conf.high = numeric(),
#                          statistic = numeric(),
#                          PValue = numeric(),
#                          stringsAsFactors = FALSE)

# Capture output to the file
capture.output({
  # Open the file connection
  sink(file_path2, append = FALSE)
  
  for (response_variable in response_variables) {
    # Fit the linear mixed-effects model with nlme
    formula_str <- paste(response_variable, "~ VFG010906", sep = "")
    formula <- as.formula(formula_str)
    
    tryCatch({
      model <- lme(formula,
                   random = list(Block = pdDiag(~1), date = pdIdent(~date_numeric)),
                   data = filtered_table_all_clr)
      
      # Print results for the current response variable
      cat("Response Variable:", response_variable, "\n")
      
      # ANOVA
      cat("ANOVA:\n")
      anova_result <- anova(model)
      df_anova_result <- data.frame(anova_result)
      df_anova_result["VFG010906","p.value"]
      
      
      # Extract ANOVA p-value and save it to the file and data frame
      anova_p_value <- anova_result$"p-value"[2]
      cat("ANOVA p-value: ", anova_p_value, "\n")
      anova_df <- rbind(anova_df, data.frame(Response_Variable = response_variable,
                                             ANOVA_PValue = anova_p_value,
                                             stringsAsFactors = FALSE))
      
      # Pairwise comparisons
      #     cat("Pairwise Comparison Results:\n")
      #     pwc <- emmeans(model, specs = "TETM", data = filtered_table_all_clr)
      #     pairwise_results <- tidy(pwc, conf.int = TRUE)
      #     print(pairwise_results)
      
      # Save pairwise results to data frame
      #     pairwise_df <- rbind(pairwise_df, data.frame(Response_Variable = response_variable,
      #                                                  Comparison = pairwise_results[, 1],
      #                                                  Estimate = pairwise_results$estimate,
      #                                                  SE = pairwise_results$std.error,
      #                                                  df = pairwise_results$df,
      #                                                  conf.low = pairwise_results$conf.low,
      #                                                  conf.high = pairwise_results$conf.high,
      #                                                  statistic = pairwise_results$statistic,
      #                                                 PValue = p.adjust(pairwise_results$p.value, method = "bonferroni"),
      #                                                 stringsAsFactors = FALSE))
      cat("\n")
    }, error = function(e) {
      # Print or log the error message
      cat(paste("Error for", response_variable, ":", conditionMessage(e), "\n"))
    })
  }
  # Close the file connection
  sink()
}, file = file_path2, append = TRUE)
# Check the file contents
cat(readLines(file_path2), sep = "\n")

# Print data frames
print(anova_df)
#print(pairwise_df)


write.csv(anova_df, file = "lme_ANOVA_pvalues_VFG010906.csv", row.names = FALSE)
#write.csv(pairwise_df, file = "lme_PWC_pvalues_VFG010906.csv", row.names = FALSE)

####Combine the two tables
anova <- read.csv("lme_ANOVA_pvalues_VFG010906.csv", header = TRUE, stringsAsFactors = FALSE) 
#pairwise <- read.csv("lme_PWC_pvalues_VFG010906.csv", header = TRUE, stringsAsFactors = FALSE) 
#merged_table_pvalues <- merge(anova, pairwise, by.x = 'Response_Variable', by.y = 'Response_Variable')
#write.csv(merged_table_pvalues, file = "lme_ALL_pvalues_VFG010906.csv", row.names = FALSE)

library(dplyr)

# Assuming you have your merged_table data frame
filtered_table <- anova %>%
  filter(ANOVA_PValue < 0.05)

# Display the result
print(filtered_table)
write.csv(filtered_table, file = "lme_ALL_pvalues_VFG010906_only_significant_ANOVA.csv", row.names = FALSE)



#############################PLS-DA##############################################
#install.packages("caret")
library(caret)
######################### create the table for PLS-DA####################################
# Assuming 'your_data' is your data frame
# and 'sites' is the response variable, and other columns are predictor variables
response_variable <- as.factor(filtered_table_all_clr$Block)
columns_to_exclude <- c("Total_per_Samples", "Total_AMR_per_Sample", "SAMPLE_ID", "Block", "date", "desc_", "Land_use", "sample_type", "Season", "sampleID", "site_yr", "yeardayofyear", "Sample_Site", "Site", "Year", "Week", "Month", "site", "strahler", "lon", "lat", "data_numeric")
# Remove specified columns
predictor_variables <- filtered_table_all_clr[, -which(names(filtered_table_all_clr) %in% columns_to_exclude)]

dim(response_variable)
dim(predictor_variables)
# Create a data frame with both the predictor and response variables
data_for_plsda <- data.frame(response_variable, predictor_variables)
class(data_for_plsda$response_variable)
data_for_plsda$response_variable <- as.factor(data_for_plsda$response_variable)
nlevels(data_for_plsda$response_variable)
# Print column names in data_for_plsda
print(colnames(data_for_plsda))
# Print AMR_VF_column_names
print(AMR_VF_column_names)
cols_exist <- AMR_VF_column_names %in% colnames(data_for_plsda)
print(cols_exist)
AMR_VF_column_names <- setdiff(AMR_VF_column_names, "SAMPLE_ID")
AMR_VF_column_names <- setdiff(AMR_VF_column_names, "Total_per_Samples")
AMR_VF_column_names <- setdiff(AMR_VF_column_names, "Total_AMR_per_Sample")
data_for_plsda[AMR_VF_column_names] <- lapply(data_for_plsda[AMR_VF_column_names], factor)
# Check the levels of each factor variable in your data
factor_levels <- sapply(data_for_plsda, function(x) if (is.factor(x)) levels(x) else NULL)
# Print only factor variables with their levels
factor_levels <- factor_levels[!sapply(factor_levels, is.null)]
print(factor_levels)
# Identify factor variables with only one level
single_level_factors <- names(factor_levels)[sapply(factor_levels, function(x) length(x) == 1)]
# Remove columns with single-level factors from the dataframe
data_for_plsda <- data_for_plsda[, !(names(data_for_plsda) %in% single_level_factors)]
# Define the control parameters for the train function
library(caret)
nzv_cols <- nearZeroVar(data_for_plsda)
data_for_plsda <- data_for_plsda[, -nzv_cols]

#install.packages("pls")
library(pls)

data_for_plsda$response_variable <- as.factor(data_for_plsda$response_variable)
str(data_for_plsda$response_variable)
summary(plsda_model)
response_matrix <- model.matrix(~ response_variable - 1, data = data_for_plsda)

# Fit PLS model
pls_model <- plsr(response_matrix ~ ., data = data_for_plsda, scale = TRUE)



# Access results
summary(pls_model)

#################################visualizations###################################################

library(ggplot2)

# Extract scores
scores <- predict(pls_model, type = "scores")

# Plot scores labelled with land use

png("Scores_Plot_Land_Use.png", width = 800, height = 600)
plot(scores[, 1], scores[, 2], col = as.numeric(as.factor(filtered_table_all_clr$Land_use)), 
     pch = 19, xlab = "Component 1", ylab = "Component 2", main = "Scores Plot of PLS-DA")
legend("topright", legend = levels(as.factor(filtered_table_all_clr$Land_use)), 
       col = 1:length(unique(filtered_table_all_clr$Land_use)), pch = 19, title = "Site")
dev.off()



# # Plot scores labeled Sample Sites
# plot(scores[, 1], scores[, 2], col = as.numeric(as.factor(filtered_table_all_clr$Block)), 
#      pch = 19, xlab = "Component 1", ylab = "Component 2", main = "Scores Plot of PLS-DA")
# legend("topright", legend = levels(as.factor(filtered_table_all_clr$Block)), 
#        col = 1:length(unique(filtered_table_all_clr$Block)), pch = 19, title = "Site")


# Plot scores labeled Season

png("Scores_Plot_Season.png", width = 800, height = 600)
plot(scores[, 1], scores[, 2], col = as.numeric(as.factor(filtered_table_all_clr$Season)), 
     pch = 19, xlab = "Component 1", ylab = "Component 2", main = "Scores Plot of PLS-DA")
legend("topright", legend = levels(as.factor(filtered_table_all_clr$Season)), 
       col = 1:length(unique(filtered_table_all_clr$Season)), pch = 19, title = "Site")
dev.off()



loadings <- pls_model$loadings[, 1:2]

# Plot loadings
png("Loading_Plot.png", width = 800, height = 600)
plot(loadings[, 1], loadings[, 2], col = 1:length(AMR_VF_column_names),
     pch = 19, xlab = "Loading Component 1", ylab = "Loading Component 2",
     main = "Loadings Plot of PLS-DA")
text(loadings[, 1], loadings[, 2], labels = AMR_VF_column_names, pos = 3, cex = 0.8)
dev.off()


############################ get all relevant information from model and save to txt file ##############################################################
# Number of components
num_components <- 240  # Change this based on the number of components in your model
output_file <- "pls_results.txt"
cat("Component\tVariable\tLoading\tWeight\n", file = output_file, append = TRUE)
# Open a file for writing
for (i in 1:num_components) {
  # Get variable loadings and weights
  loadings <- pls_model$loadings[, i]
  weights <- pls_model$weights[, i]
  
  
  # Print header for the current component
  cat(paste("Component", i, "\n", sep = "\t"), file = output_file, append = TRUE)
  
  # Print results for the current component
  cat(paste("Variable\tLoading\tWeight\n", sep = "\t"), file = output_file, append = TRUE)
  for (j in seq_along(loadings)) {
    cat(paste(names(loadings)[j], loadings[j], weights[j, 1], sep = "\t"), "\n", file = output_file, append = TRUE)
  }
  
  # Add an empty line between components
  cat("\n", file = output_file, append = TRUE)
}


#Just getting variance explained by each component
# Extract explained variance
options(max.print = 1000000)
summary_output <- capture.output(summary(pls_model))

# Open a file for writing
output_file <- "for_samples_components_explained_variance.txt"
writeLines(summary_output, output_file)


