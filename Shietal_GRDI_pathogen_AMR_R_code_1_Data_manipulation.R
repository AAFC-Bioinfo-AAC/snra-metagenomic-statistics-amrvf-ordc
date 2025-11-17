############################################################################
#####################Amplicon data--whole bacterial community #################
##########################################################################


###loading feature table and taxonomy table from amplicon data pipeline
asv.dada2 <- read.table("feature-table.tsv", header = TRUE, sep = "\t")
taxa.dada2 <- read.table("taxonomy.tsv", header = TRUE, sep = "\t")

#merge to asv table
dim(asv.dada2) #[1] 307975   1519
names(asv.dada2) #OTU.ID
names(taxa.dada2)[1] <- "OTU.ID"
names(taxa.dada2)[2] <- "taxonomy"
rownames(asv.dada2) <- asv.dada2$OTU.ID
asv.dada2 <- asv.dada2[,-1]
asv_table <- merge(asv.dada2, taxa.dada2, by = "OTU.ID")
dim(asv_table)

##meta data
sample.list <- read.csv("sample.list_16_21.csv", header = TRUE, sep =",")
identical(names(asv.dada2),sample.list$sampleID)
dim(asv.dada2)
dim(sample.list)



#Exclude the samples with total reads less than 1000

asv.dada2.f <- asv.dada2[,colSums(asv.dada2) >= 1000]
dim(asv.dada2.f)
rownames(taxa.dada2) <- taxa.dada2$OTU.ID
asv.dada2.f$OTU.ID <- rownames(asv.dada2.f)

asv_table <- left_join(asv.dada2.f, taxa.dada2, by = "OTU.ID")

#sample.list$sampleID <- names(asv.dada2)


#the current asv table include samples filtered by 0.2 um and 0.7
#Now need subset the data to only include Filter == 0.2


# Filter the sample list for Filter == 0.2
sample.0.2 <- sample.list[sample.list$Filter == 0.2,]
sample.0.2 <- sample.0.2[1:693,]
filtered_samples <- sample.0.2$sampleID

sample.soil <- sample.list[!sample.list$Source == "Surfwat",]
filtered_samples_2 <- sample.soil$sampleID

asv.dada2.0.2 <- asv.dada2.f[, names(asv.dada2.f) %in% filtered_samples]
asv.dada2.soil <- asv.dada2[,names(asv.dada2) %in% filtered_samples_2]

# Check the resulting subset
head(asv.dada2.0.2)
dim(asv.dada2.0.2) #307975    680
names(asv.dada2.0.2)
rownames(asv.dada2.0.2)
asv.dada2.0.2$OTU.ID <- rownames(asv.dada2.0.2)
asv_table.0.2 <- merge(asv.dada2.0.2, taxa.dada2, by = "OTU.ID")


asv.dada2.soil$OTU.ID <- rownames(asv.dada2.soil)
asv_table.soil <- merge(asv.dada2.soil, taxa.dada2, by = "OTU.ID")
min(asv_table.0.2$Confidence)
max(asv_table.0.2$Confidence)

asv_table_0.2_filter <- asv_table.0.2 %>%
  filter(Confidence > 0.7) #293145    683

asv_table_0.2_filter <- asv_table_0.2_filter  %>%
  filter(rowSums(asv_table_0.2_filter[,2:681])>0) #166747    683
rownames(asv_table_0.2_filter) <- asv_table_0.2_filter$OTU.ID

asv_table_soil_filter <- asv_table.soil %>%
  filter(Confidence > 0.7) 

asv_table_soil_filter <- asv_table_soil_filter  %>%
  filter(rowSums(asv_table_soil_filter[,2:370])>0) #95492   372


###Change the column name with the format as "SN024.20190624.Surfwat", check if there are repeated names, combine them with sum.

# Step 1: Rename the columns as before
id_mapping <- setNames(sample.list$ID_combined, sample.list$sampleID)
colnames_to_modify <- colnames(asv_table_0.2_filter)[2:681]
new_colnames <- ifelse(colnames_to_modify %in% names(id_mapping), 
                       id_mapping[colnames_to_modify], 
                       colnames_to_modify)
colnames(asv_table_0.2_filter)[2:681] <- new_colnames

asv_table_0.2_filter <- asv_table_0.2_filter[,-ncol(asv_table_0.2_filter)] #166747    683


# Step 2: Handle duplicates by summing rows
# S1: Extract the first 22 characters of column names to define groups
col_groups <- substr(colnames(asv_table_0.2_filter), 1, 22)

# Step 2: Create a data frame to map col_groups to original column names
col_map <- data.frame(
  OriginalName = colnames(asv_table_0.2_filter),
  GroupName = col_groups,
  stringsAsFactors = FALSE
)

# Step 3: Transpose the data frame excluding metadata (e.g., OTU.ID or taxonomy)
metadata_cols <- c("OTU.ID","taxonomy")  # Adjust if different metadata column names exist
data_cols <- setdiff(colnames(asv_table_0.2_filter), metadata_cols)
asv_transposed <- as.data.frame(t(asv_table_0.2_filter[, data_cols]))

# Step 4: Aggregate columns by summing them using their GroupName
aggregated <- rowsum(asv_transposed, group = col_groups[-c(1,682)])

# Step 5: Transpose back to original orientation
asv_summed <- t(aggregated)

# Step 6: Reattach metadata columns (e.g., OTU.ID and taxonomy)
asv_table_0.2_filter_final <- cbind( asv_summed,asv_table_0.2_filter[, metadata_cols, drop = FALSE]) #166747    563

rownames(asv_table_0.2_filter_final) 

# Check the final data frame
dim(asv_table_0.2_filter_final)




####Filter the samples to correspond with mm.499
names(asv_table_0.2_filter_final) <- gsub(".Surfwat","",names(asv_table_0.2_filter_final))
asv_table_16S_allsample <- asv_table_0.2_filter_final[,-562]

asv_table_16S_mm499 <- asv_table_16S_allsample[,names(asv_table_16S_allsample) %in% rownames(mm.499)]
asv_table_16S_mm499 <- cbind(asv_table_16S_mm499, asv_table_0.2_filter_final$taxonomy)
names(asv_table_16S_mm499)[500] <- "taxonomy"

asv_table_16S_sgsample <- asv_table_0.2_filter_final[,names(asv_table_0.2_filter_final) %in% names(mm.499)]
asv_table_16S_sgsample <- asv_table_16S_sgsample[rowSums(asv_table_16S_sgsample[,-ncol(asv_table_16S_sgsample)])>0,] #107175    271
dim(asv_table_16S_sgsample)
asv_table_16S_sgsample_g <- RAM::tax.abund(asv_table_16S_sgsample, rank = "g")


xx <- names(asv.table.shotgun[!(names(asv.table.shotgun) %in% names(asv_table_16S_sgsample))]) 
xx
#"SN006.20201102" "SN018.20201102" "SN018.20210628" "SN019.20170704" "SN019.20201102" "SN020.20201019" "SN021.20201102" "SN024.20160912" "SN024.20201102"
list1 <- c("SN006.20201102", "SN018.20201102", "SN019.20170704", "SN019.20201102", "SN020.20201019", "SN021.20170717", "SN021.20180709", "SN021.20201102", "SN024.20160912", "SN024.20201102")

rownames(mm1)
mm.f <- mm1  %>% 
  filter(!rownames(mm1) %in% list1) #272 176
dim(mm1)


##########################################
#########Create Phyloseq object#########
##########################################

##Creat phyloseq file
asv <- asv_table_16S_allsample[,-ncol(asv_table_16S_allsample)] # change to asv_table_16S_sgsample for 270 samples
asv <- asv[,names(asv) %in% rownames(mm)] #166745    499
asv_table_16S_allsample <- cbind(asv, asv_table_16S_allsample$taxonomy) #166745    500
asv_table_16S_allsample <- asv_table_16S_allsample[rowSums(asv_table_16S_allsample[,-ncol(asv_table_16S_allsample)])>0,] #152076    500
names(asv_table_16S_allsample)[500] <- "taxonomy"
tax <- asv_table_16S_allsample %>%
  separate(
    taxonomy, 
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
    sep = "; ", 
    remove = FALSE
  ) 

mm <- mm[rownames(mm) %in% names(asv),]
head(tax)
# Extract the taxonomy table
tax <- tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species","taxonomy")]
tax <- tax %>%
  mutate(across(1:7, ~ sub("^[a-z]__", "", .)))

head(tax)

asv_mat <- as.matrix(asv)
taxa_mat <- as.matrix(tax)
OTU = otu_table(asv_mat, taxa_are_rows = TRUE)
TAX = tax_table(taxa_mat)
samples = sample_data(mm)

##whole community of 16S data for 499 samples
physeq.16s.499 <-  phyloseq(OTU, TAX, samples)

##whole community of 16S data for 270 samples, corresponding to shotgun data samples
physeq.16s.270 <-  phyloseq(OTU, TAX, samples)







#####Generate genus table from phyloseq######
# Step 1: Aggregate to genus level
phyloseq_genus <- tax_glom(physeq.16s.270, taxrank = "Genus")

# Step 2: Extract OTU table (genus-level abundance table)
genus_table <- as.data.frame(otu_table(phyloseq_genus))

# Step 3: Ensure samples are rows and genera are columns
if (taxa_are_rows(phyloseq_genus)) {
  genus_table <- t(genus_table)  # Transpose if needed
}

# Step 4: Add genus names as column names (optional)
tax_table_genus <- as.data.frame(tax_table(phyloseq_genus))
colnames(genus_table) <- tax_table_genus$Genus
genus_table <- as.data.frame(genus_table)
# Step 5: Output genus table
head(genus_table)  # Display first rows






########################################################################
#####################Shotgun data--whole bacterial community #################
######################################################################


#read asv and tax information

library(data.table)
s.o <- readr::read_tsv("species_matrix.tsv")
s.o <- as.data.frame(s.o)

genus <- readr::read_csv("genus_cross_table.csv", col_names = TRUE)
genus <- genus %>% mutate(across(everything(), as.numeric))

genus <- t(genus)
genus <- as.data.frame(genus)
colnames(genus) <- genus[1,]
genus <- genus[-1,]
names(genus) <- gsub("-",".",names(genus))
genus <- as.data.frame(genus)
genus.270 <- genus[,names(genus) %in% row.names(genus_table)]

genus.270 <- t(genus.270)
genus.270 <- as.data.frame(genus.270)
#To convert all columns of a data frame to numeric
library(dplyr)

genus.270 <- genus.270 %>% mutate(across(everything(), as.numeric))


genus.270 <- genus.270[,colSums(genus.270)>0]
matching_columns_sg <- names(genus.270)[str_detect(names(genus.270), str_c(fecal_combined_genus, collapse = "|"))]

genus.fecal <- genus.270[,matching_columns_sg]
genus.fecal[] <-  lapply(genus.fecal, function(x) as.numeric(as.character(x)))
str(genus.fecal)
genus.fecal <- genus.fecal[,colSums(genus.fecal)>0]
dim(genus.fecal)
rowSums(genus.fecal)

genus.fecal.ra <- decostand(genus.fecal, "total",MARGIN = 1)
genus.sg.pathogen <-genus.fecal[,names(genus.fecal) %in% pathogen.list]
dim(genus.sg.pathogen)
genus.sg.pathogen.ra <- genus.fecal.ra[,names(genus.fecal.ra) %in% pathogen.list]

s.o$taxonomy <- s.o$TAXONOMY

s.o <- s.o %>%
  separate(taxonomy, 
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "\\|",
           remove = FALSE) %>%  # Keep original column if desired
  mutate(across(everything(), ~ sub("^[a-z]__", "", .)))  # Remove prefixes (e.g., k__, p__)
#names(s.o)[328] <- "taxonomy"
# We find repeat samples with sample names of SN020-20201019 & SN20-20201019, and SN020-20201102&SN20-20201102
#I decide to delete the two samples SN20-20201019 and SN20-20201102
s.o <- s.o[, !(colnames(s.o) %in% c("SN20-20201019", "SN20-20201102"))]
s.o <- as.data.frame(s.o)
rownames(s.o) <- paste("asv", seq_len(nrow(s.o)), sep = "")

# Extract the taxonomy table
tax <- s.o[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species","TAXONOMY")]
rownames(tax)

# Extract the asv table
asv <- s.o[,-c(326:333)]

asv <- asv[,-1]
rownames(asv)
names(asv) <- gsub("-",".",names(asv))

mm.div.16S.fecal$sampleID <- gsub(".Surfwat","",mm.div.16S.fecal$sampleID)
# Only keep the samples in asv we found in the meta data file
asv <- asv[,names(asv) %in% mm.div.16S.fecal$sampleID]
dim(asv)
dim(mm1)
identical(names(asv), rownames(mm.div.16S.fecal))

#change all of columns of asv to numeric
asv[] <- lapply(asv, function(x) as.numeric(as.character(x)))

#remove the asv in asv file that absent across all the samples after removing some samples
asv <- asv[rowSums(asv)>0, ]
dim(asv)
rownames(asv)
#remove the asv from tax file
tax <- tax[rownames(tax) %in% rownames(asv),]

mm.270 <- mm.div.16S.fecal[,1:176]

##Creat phyloseq file
asv_mat <- as.matrix(asv)
ta_mat <- as.matrix(tax)
OTU = otu_table(asv_mat, taxa_are_rows = TRUE)
TAX = tax_table(ta_mat)
samples = sample_data(mm.270)

physeq.shotgun.270 <-  phyloseq(OTU, TAX, samples)

saveRDS(physeq.shotgun.270, file = "physeq.shotgun.270.rds")


########################################################
###########Fecal-associated genera list#############
#####################################################


#########Step 1--fecal-community list compling#######
amdb <- read.csv("AMDB.csv", header = TRUE, sep = ",")
dim(amdb)
addagma <- read.csv("ADDAGMA_genus.csv", header = TRUE, sep = ",")
amdb_genus <- amdb$Genus
addagma_genus <- addagma$Genus

fecal_combined_genus <- unique(c(amdb_genus, addagma_genus)) #combine two datasets of genera list: amdb and addagma
write.csv(fecal_combined_genus,"fecal_combined_genus.csv")



###########Step2: Fecal-associated community--16S#############
#

# Check that your taxonomic rank is "Genus"
# Replace "Genus" with your actual taxonomic rank name if different
taxa_genus <- tax_table(physeq.16s.499)[, "Genus"]
taxa_genus <- as.data.frame(taxa_genus)
g.list <- unique(taxa_genus$Genus)
matching_columns_sg <- g.list[g.list %in% fecal_combined_genus]

# Find taxa that match your list
matching_taxa <- taxa_names(physeq.16s.499)[taxa_genus[,1] %in% matching_columns_sg]

# Subset the phyloseq object to only include those taxa
physeq.16S.fecal.499 <- prune_taxa(matching_taxa, physeq.16s.499)

physeq.16S.fecal.270 <- prune_samples(rownames(mm.270), physeq.16S.fecal.499)
sample_data(physeq.16S.fecal.270) <- sample_data(mm.270)
physeq.16S.fecal.270
physeq.16S.fecal.270 <- prune_taxa(taxa_sums(physeq.16S.fecal.270) > 0, physeq.16S.fecal.270)

saveRDS(physeq.16S.fecal.270, file = "physeq.16S.fecal.270.rds")

saveRDS(physeq.16S.fecal.499, file = "physeq.16S.fecal.499.rds")


###########Step 3: Fecal-associated community--shotgun#############


# Check that your taxonomic rank is "Genus"
# Replace "Genus" with your actual taxonomic rank name if different
taxa_genus <- tax_table(physeq.shotgun.270)[, "Genus"]
taxa_genus <- as.data.frame(taxa_genus)
g.list <- unique(taxa_genus$Genus)
matching_columns_sg <- g.list[g.list %in% fecal_combined_genus]

# Find taxa that match your list
matching_taxa <- taxa_names(physeq.shotgun.270)[taxa_genus[,1] %in% matching_columns_sg]

# Subset the phyloseq object to only include those taxa
physeq.shotgun.fecal.270 <- prune_taxa(matching_taxa, physeq.shotgun.270)

saveRDS(physeq.shotgun.fecal.270, file = "physeq.shotgun.fecal.270.rds")


#################################################
#######Pathogen community_16S##################
#############################################

#######Step 1--Pathogen list compling#########




###Zoonotic or opportunisitic pathogens
#ref1: Compiled pathogen list infecting humans (Bartlett et al., 2022 )
pathogen.human <- read.csv("Pathogen_Bartlett_2022_full list.csv",  header = TRUE, sep = ",")

#ref2:NCBI Pathogen Detection database (https://www.ncbi.nlm.nih.gov/pathogens/)
pathogen.ncbi <- read.csv("NCBI_pathogen_detection_list.csv", header = TRUE, sep = ",")

#ref3:Voluntary 2025 U.S. National Animal Health Reporting System (NAHRS) Reportable Diseases, Infections, and Infestations List 
pathogen.nahrs <- read.csv("NHARS_bacteria_disease_list.csv", header = TRUE, sep = ",")

#ref4: PATRIC pathogen database (https://www.bv-brc.org/)

PATRIC_list <- c("Acinetobacter", "Bacillus",
                 "Bartonella", "Borreliella","Brucella",
                 "Burkholderia", "Campylobacter", "Chlamydia","Clostridium","Coxiella",
                 "Ehrlichia","Escherichia", "Francisella", "Helicobacter","Listeria", "Mycobacterium",
                 "Pseudomonas","Rickettsia","Salmonella","Shigella","Staphylococcus","Streptococcus",
                 "Vibrio","Yersinia", "Neisseria", "Enterococcus")

#ref4: CDC zoonotic pathogens: https://www.cdc.gov/healthy-pets/diseases/index.html

pathogen.cdc <- read.csv("CDC_bacterial pathogen list.csv", header = TRUE, sep = ",")
CDC <- unique(pathogen.cdc$Genus)

#ref5: WHO bacterial priority pathogens list, 2024 
#(https://www.who.int/publications/i/item/9789240093461)


pathogen.who <- read.csv("WHO_bacteria_list.csv", header = TRUE, sep = ",")
WHO <- pathogen.who$Genus

xx <- c(CDC, WHO)


pathogen.list <- unique(c(pathogen.human$genus, pathogen.ncbi$Genus, pathogen.nahrs$Genus, PATRIC_list,CDC, WHO))

write.csv(pathogen.list, "output_pathogen.list.genera.csv")



#######Step 2--Pathogen community--16S#########
physeq.16S.fecal.270

physeq.genus <- tax_glom(physeq.16S.fecal.270, taxrank = "Genus")
genus.16S.fecal.mm270 <- otu_table(physeq.genus)
class(genus.16S.fecal.mm270)
genus.16S.fecal.mm270 <- t(genus.16S.fecal.mm270)
dim(genus.16S.fecal.mm270)

genus.16S.fecal.mm270 <- as.data.frame(as.matrix(genus.16S.fecal.mm270))
genus.names <- as.character(tax_table(physeq.genus)[, "Genus"])
colnames(genus.16S.fecal.mm270) <- genus.names




library("stringr")
matching_columns4 <- names(genus.16S.fecal.mm270)[str_detect(names(genus.16S.fecal.mm270), str_c(pathogen.list, collapse = "|"))]
genus.16S.pathogen.270 <- genus.16S.fecal.mm270[, matching_columns4]

tax <- tax_table(physeq.16S.mm270)
colnames(tax)

pathogen.taxa <- taxa_names(physeq.16S.mm270)[tax[, "Genus"] %in% matching_columns4]
physeq.16S.pathogen.270 <- prune_taxa(pathogen.taxa, physeq.16S.mm270)
physeq.16S.pathogen.270 <- prune_taxa(taxa_sums(physeq.16S.pathogen.270) > 0, physeq.16S.pathogen.270)
physeq.16S.pathogen.270

saveRDS(physeq.16S.pathogen.270, file = "physeq.16S.pathogen.270.rds")



#######Step 3--Pathogen community--shotgun#########
physeq.shotgun.fecal.270 <- readRDS("physeq.shotgun.fecal.270.rds")


# Extract tax_table as a data.frame
tax_df <- as.data.frame(tax_table(physeq.shotgun.fecal.270))

# Check if Genus column exists (adjust "Genus" if your column has a different name)
table(tax_df$Genus %in% pathogen.list)

# Identify taxa indices matching pathogen genera
taxa_to_keep <- rownames(tax_df)[tax_df$Genus %in% pathogen.list]

# Prune taxa in the phyloseq object
physeq.shotgun.pathogen.270 <- prune_taxa(taxa_to_keep, physeq.shotgun.fecal.270)
saveRDS(physeq.shotgun.pathogen.270, "physeq.shotgun.pathogen.270.rds")



library("stringr")
matching_columns4 <- names(genus.16S.fecal.mm270)[str_detect(names(genus.16S.fecal.mm270), str_c(pathogen.list, collapse = "|"))]
genus.16S.pathogen.270 <- genus.16S.fecal.mm270[, matching_columns4]

tax <- tax_table(physeq.16S.mm270)
colnames(tax)

pathogen.taxa <- taxa_names(physeq.16S.mm270)[tax[, "Genus"] %in% matching_columns4]
physeq.16S.pathogen.270 <- prune_taxa(pathogen.taxa, physeq.16S.mm270)
physeq.16S.pathogen.270 <- prune_taxa(taxa_sums(physeq.16S.pathogen.270) > 0, physeq.16S.pathogen.270)
physeq.16S.pathogen.270


saveRDS(physeq.16S.pathogen.270, file = "physeq.16S.pathogen.270.rds")
