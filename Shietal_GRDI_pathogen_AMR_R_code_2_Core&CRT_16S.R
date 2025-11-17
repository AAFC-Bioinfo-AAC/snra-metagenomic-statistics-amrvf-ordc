#############################################
###########1-Data preparation############
##############################################

## Part 1: Identification CRT##
#####identification of CRT for 499 samples across all sites 
#Prepare dataset

library(vegan)
physeq.16S.fecal.499 <- readRDS("physeq.16S.fecal.499.rds")

asv.16S.all.fecal <- as.data.frame(otu_table(physeq.16S.fecal.499))
dim(asv.16S.all.fecal)

# Set a rarefaction depth (choose based on min sample size)

asv <- asv.16S.all.fecal[,colSums(asv.16S.all.fecal) > 1000]
dim(asv)

# Perform rarefaction
min_depth <- min(colSums(asv)) 
asv.rf <- rrarefy(t(asv), sample = min_depth)
asv.rf <- t(asv.rf)
asv.rf <- as.data.frame(asv.rf[rowSums(asv.rf)>0,])
mm <- mm.499[rownames(mm.499) %in% names(asv.rf),] #487 samples
mm$sampleID <- gsub(".Surfwat","",mm$sampleID)
mm$YW <- paste0(mm$Year,"_",mm$Week)
mm$YW <- as.factor(mm$YW)
mm.a <- mm %>%
  filter(site_type == "agr")
mm.m <- mm %>%
  filter(site_type == "mixed")
mm.f <- mm %>%
  filter(site_type == "forest")
mm.ac <- mm.a %>%
  filter(site %in% c("18", "19"))
mm.am <- mm.a %>%
  filter(site %in% c("20", "21"))

# Ensure your OTU table is a dataframe or tibble
otu.a <- asv.rf[,names(asv.rf) %in% rownames(mm.a)]
otu.a <- otu.a[rowSums(otu.a)>0,]
otu.m <- asv.rf[,names(asv.rf) %in% rownames(mm.m)]
otu.m <- otu.m[rowSums(otu.m)>0,]

otu.ac <- asv.rf[,names(asv.rf) %in% rownames(mm.ac)]
otu.ac <- otu.ac[rowSums(otu.ac)>0,] #4894

otu.am <- asv.rf[,names(asv.rf) %in% rownames(mm.am)]
otu.am <- otu.am[rowSums(otu.am)>0,] #4429

write.table(asv.rf, "16s_asv_rf.txt", sep = "\t", quote = FALSE)
write.table(otu.a, "asv_rf_agr.txt", sep = "\t", quote = FALSE)
write.table(otu.m, "asv_rf_mixed.txt", sep = "\t", quote = FALSE)

#############################################
###########2-Conditionally rare taxa############
##############################################

##https://github.com/ShadeLab/ConditionallyRareTaxa
#Step 1.
#If they are not installed already, install the following required R packages: vegan, TSA.  Then, load the libraries to the R workspace by copying and pasting the commands below into the R console:
#library(vegan)
library(TSA)

#Step 2.
#Place the input file and script in the same working directory to run this script.  Change the working directory in R to match where the files have been placed.

#Step 3.
#Load the necessary functions into your R workspace, contained in a separate file, "CRT_functions.R" 
source("CRT_Functions_v1.1.R")

#Step 4.  
#Change the options below to match your dataset.  The options are:  
#otu_fp - type the the full name of your dataset file, including the extension
#abund_thresh -  Change the maximum abundance threshold, if desired. Defaults to 0.005
#abund_thresh_ALL - Use TRUE if you want to use the full dataset (ALL OTUs) to calculate relative abundances.  Use FALSE if you want to use the non-singleton (filtered) dataset to calculate relative abundances.  Default is FALSE.
#b_thresh - Change the coefficient of bimodality threshold, if desired.  Defaults to 0.90
#rdp_lastcol - Use TRUE if the last column of the dataset contains the taxonomic assignments of OTUs, use FALSE if not
#Then,to run the script, copy and paste the command into the R console:

###We need remove space 
SimpleRareToPrev.f(otu_fp="16s_asv_rf.txt",abund_thresh=0.005, abund_thresh_ALL=FALSE,b_thresh=0.90, rdp_lastcol=TRUE)
#we could run asv_rf.txt for all sites,
#asv_rf_agr for agri_ditches sites and 
#asv_rf_mixed for mixed sites

####obtain the asv names identified for CRT 
##remove the space in taxonomy column of the result

crt.all <- read.table("ResultsFile_ConditionallyRareOTUID_0.005_0.9_NOSIG.txt", header = TRUE, sep = "\t") # remove the space between ranks in the file.
names(crt.all)
crt.list <- crt.all$OTUID



asv.crt <- asv.rf[rownames(asv.rf) %in% crt.list,] 
dim(asv.crt) #1385
#replace otu by otu.a.d for specific CRT for ad


asv.crt.n <- asv.crt[!(rownames(asv.crt) %in% rownames(asv.core)),] 
dim(asv.crt.n) #1383


write.csv(asv.crt.n, file="asv.crt.n_new.csv", quote = FALSE)

asv.crt.n <- read.csv("asv.crt.n_new.csv", header= TRUE, sep = ",")
rownames(asv.crt.n) <- asv.crt.n$ASV_ID
asv.crt.n <- asv.crt.n[,-1]

#####################################################
###### 3-Identification of core microbiome#######################
######################################################

library(tidyverse)
library(reshape2)
##Prioritizing core microbiome based on TIME
#######

######Calculate Index
#a.d land use types
nReads=1018      # input dataset needs to be rarified and the rarifaction depth included 


levels(factor(mm$YW))
counts <- table(factor(mm$YW)) #count sample numbers in each level

# Step 2: Keep levels with more than one sample
valid_levels <- names(counts[counts > 1])#58 for agr; 71 for mixed;73 for all samples

asv<- asv.rf
col <- 73 # number of time points levels
row <- 11146 #number of asv
mm.x <- mm

#Create a matrix for recording results of time-specific occupancy
OTU_acc <-array(-99,dim=c(row,col)) #asv numbers and levels of time points exclude the levels only one replicate
rownames(OTU_acc) <- rownames(asv)
colnames(OTU_acc) <- valid_levels
#a <- as.factor(mm.rf$Date)
#count(mm.rf$Date)
dim(OTU_acc)
##Create a matrix for recording results of replication consistency
OTU_rep <-array(-99,dim=c(row,col))
rownames(OTU_rep) <- rownames(asv)
colnames(OTU_rep) <- valid_levels
dim(OTU_rep)
##Create a matrix for recording results of index
OTU_index <- array(-99,dim=c(row,3))
rownames(OTU_index) <-rownames(asv)
colnames(OTU_index) <- c("sumF", "sumG", "index")
dim(OTU_index)

for ( date in valid_levels)  {
  names <-rownames(mm.x)[as.character(mm.x$YW) %in% date]
  asv.sub <- asv[,names]
  asv.sub <- asv.sub[rowSums(asv.sub)>0,]
    for (o in rownames(asv.sub)) {
      asv.sub_o <- asv.sub[o,]
      asv.sub_o_PA <- decostand(asv.sub_o, method = "pa") 
    
    time_freq=sum(asv.sub_o_PA)/ncol(asv.sub_o_PA) # frequency of detection between time points
    coreTime=ifelse(time_freq == 1, 1, 0) # 1 only if occupancy 1 with specific time, 0 if not
    
    OTU_acc[o,date]=time_freq
    OTU_rep[o,date]=coreTime
    
  }
}


OTU_acc[OTU_acc < 0] <- 0
OTU_rep[OTU_rep < 0] <- 0
OTU_index[,"sumF"] <- rowSums(OTU_acc)
OTU_index[,"sumG"] <- rowSums(OTU_rep)
OTU_index[,"index"] <- (as.numeric(OTU_index[,"sumF"])+ as.numeric(OTU_index[,"sumG"]))/(ncol(OTU_acc)*2)    # calculating weighting Index based on number of time points detected
OTU_index <- as.data.frame(OTU_index)
OTU_index <- OTU_index[order(OTU_index[,3],  decreasing = TRUE),]  

#as.data.frame(sapply(OTU_index, function(x) gsub("\"", "", x)))
dim(OTU_index)
#OTU_index <- OTU_index[,-c(5:6)]
OTU_index$otu <- rownames(OTU_index)
# OTU_index$index <- as.numeric(OTU_index$index)
OTU_index$dd <- 1/(as.numeric(OTU_index[,3]))
OTU_index <- OTU_index %>%
  mutate(rank_order = dense_rank(dd)) %>%
  arrange(rank_order) 


write.csv(OTU_index, "all_index_allsamples.csv", quote=FALSE)




# Calculating the contribution of ranked OTUs to the BC similarity
BCaddition <- NULL

# calculating BC dissimilarity based on the 1st ranked OTU
otu_start=OTU_index$otu[which(OTU_index$rank_order %in% 1)]                  
start_matrix <- asv[otu_start,]
#start_matrix <- t(start_matrix[,-ncol(start_matrix)])
x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(2*nReads))
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
df_s <- data.frame(x_names,x)
names(df_s)[2] <- 1 
BCaddition <- rbind(BCaddition,df_s)
# calculating BC dissimilarity based on additon of ranked OTUs from 2nd to 500th. Can be set to the entire length of OTUs in the dataset, however it might take some time if more than 5000 OTUs are included.
for(i in 2:200){  #all:2945                         
  otu_add=OTU_index$otu[which(OTU_index$rank_order %in% i)]                        
  add_matrix <- asv[otu_add,]
  #dim9a#add_matrix <- t(add_matrix)
  start_matrix <- rbind(start_matrix, add_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_a <- data.frame(x_names,x)
  names(df_a)[2] <- i 
  BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))
}
# calculating the BC dissimilarity of the whole dataset (not needed if the second loop is already including all OTUs) 
otu <- asv[,]
x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))   
x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
df_full <- data.frame(x_names,x)
names(df_full)[2] <- length(rownames(otu))
BCfull <- left_join(BCaddition,df_full, by='x_names')

rownames(BCfull) <- BCfull$x_names
#temp_BC <- BCfull
#temp_BC$x_names <- NULL
#temp_BC_matrix <- as.matrix(temp_BC)

temp_BC <- BCaddition
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%            # mean Bray-Curtis dissimilarity
  arrange(desc(-MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))   # proportion of the dissimilarity explained by the n number of ranked OTUs
Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
BC_ranked <- left_join(BC_ranked, increaseDF)
BC_ranked <- BC_ranked[-nrow(BC_ranked),]
dim(BC_ranked)
names(BC_ranked)

write.table(BC_ranked, "BC_ranked_allsamples.txt", quote= FALSE)
#library(data.table)
#BC_ranked <- BC_ranked[,1:4]
#BC_ranked <- BC_ranked[-2578,]
write.csv(BC_ranked, file ="all_BC_ranked_allsamples.csv", quote= FALSE) #data.table


#Creating thresholds for core inclusion 
#BC_ranked$Increase_5p <- rollmean(BC_ranked$IncreaseBC, k=5, fill = NA) #package of zoo
#core  ## the last increase of 1% in the first 100 ranks; For Agriculture sites: rank = 79; for mixed site: rank = 60
#Method: 
#A) Elbow method (first order difference) (script modified from https://pommevilla.github.io/random/elbows.html)
#fo_difference <- function(pos){
#  left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
#  right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
#  return(left - right)
#}
#BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

#elbow <- which.max(BC_ranked$fo_diffs)

BC_ranked1 <- BC_ranked[1:100,]
lastCall <- last(BC_ranked1$rank[(BC_ranked1$IncreaseBC>=1.01)]) #63 (2%)

#Creating plot of Bray-Curtis similarity
ggplot(BC_ranked[1:131,], aes(x=factor(BC_ranked$rank[1:131], levels=BC_ranked$rank[1:131]))) +
  geom_point(aes(y=proportionBC)) +
  theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
  geom_vline(xintercept=last(BC_ranked1$rank[(BC_ranked1$IncreaseBC>=1.02)]), lty=3, col='blue', cex=.5) +
  labs(x='ranked OTUs',y='Bray-Curtis similarity') +
  annotate(geom="text", x=last(BC_ranked1$rank[(BC_ranked1$IncreaseBC>=1.01)]), y=.4, label=paste("Last 1% increase (",last(BC_ranked1$rank[(BC_ranked1$IncreaseBC>=1.02)]),")",sep=''), color="blue")
#rank:98 (2%);  

##based on ad samples
OTU_index <- read.csv("all_index_allsamples.csv", header = TRUE, sep = ",")
core.list <- OTU_index[OTU_index$rank_order %in% c(1:80),] 
core.list <- core.list$otu

asv.core <- asv.rf[rownames(asv.rf) %in% core.list,] 
rownames(asv.core)
dim(asv.ad.core)#135 ASVs 
dim(asv.m.core)#177 ASVs 
dim(asv.core) #98 ASVs


write.csv(asv.core, file = "asv.core_new.csv", quote = FALSE)
###


#############################
######4-top20 genera barplot
##############################

#-------------------------
#p.top20.core---------
#---------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

taxonomy_col <- "Taxonomy"

# Extract genus (adjust pattern depending on your taxonomy format)
asv.core <- asv.core %>%
  mutate(Genus = sub(".*g__", "", !!sym(taxonomy_col)))

# Extract genus from taxonomy (adjust if different format)
asv.crt.n <- asv.crt.n %>%
  mutate(Genus = sub(".*g__", "", !!sym(taxonomy_col)))


# Extract all genera across both datasets
genus.core <- unique(asv.core$Genus)
genus.crt <- unique(asv.crt.n$Genus)

# Combine all unique genera (including "Other" if present)
all.genera <- sort(unique(c(genus.core, genus.crt, "Other")))

library(viridis)

# Create color palette with enough distinct colors
high_contrast_colours <- setNames(
  viridis(length(all.genera), option = "turbo"),
  all.genera
)





# Convert to long format (ASVs as rows, samples as columns)
asv.long <- asv.core %>%
  pivot_longer(
    cols = -c(Genus, !!sym(taxonomy_col)),
    names_to = "Sample",
    values_to = "Abundance"
  )

# Add metadata (site_type)
asv.long <- asv.long %>%
  left_join(mm.499 %>% select(Sample, site_type), by = "Sample")

# Aggregate by Genus and site_type
genus.sum <- asv.long %>%
  group_by(site_type, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(site_type) %>%
  mutate(RelAbundance = Abundance / sum(Abundance))

# Optional: order site_type levels
genus.sum$site_type <- factor(genus.sum$site_type, levels = c("Forest", "Mixed", "Agriculture"))

# Plot stacked bar by site_type
p.top20.core <- ggplot(genus.sum.core, aes(x = site_type, y = RelAbundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = high_contrast_colours) +
  labs(
    x = "Site Type",
    y = "Relative Abundance (%)",
    fill = "Genus"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )



#-------------------------
#p.top20.crt---------
#---------------------------




# Convert to long format
asv.long <- asv.crt.n %>%
  pivot_longer(
    cols = -c(Genus, !!sym(taxonomy_col)),
    names_to = "sampleID",
    values_to = "Abundance"
  )

# Merge metadata info (Sample ↔ site_type)
asv.long <- asv.long %>%
  left_join(mm.499 %>% select(sampleID, site_type), by = "sampleID")

genus.sum <- asv.long %>%
  group_by(site_type, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(site_type) %>%
  mutate(RelAbundance = Abundance / sum(Abundance))

top20 <- genus.sum %>%
  group_by(Genus) %>%
  summarise(Total = sum(RelAbundance)) %>%
  top_n(20, Total) %>%
  pull(Genus)

genus.sum <- genus.sum %>%
  mutate(Genus = ifelse(Genus %in% top20, Genus, "Other")) %>%
  group_by(site_type, Genus) %>%
  summarise(RelAbundance = sum(RelAbundance), .groups = "drop")

p.top20.crt <- ggplot(genus.sum.crt, aes(x = site_type, y = RelAbundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = high_contrast_colours) +
  labs(
    x = "Site Type",
    y = "Relative Abundance (%)",
    fill = "Genus"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )



p.top20.corecrt <- plot_grid(p.top20.core + ggtitle("(A) Core"),
                             p.top20.crt + ggtitle("(B) CRT"), ncol = 2)

pdf("Supplementary Fig. S3- Top20genera_coreCRT.pdf", width = 12, height = 6)
p.top20.corecrt
dev.off()

tiff("Supplementary Fig. S3- Top20genera_coreCRT.tiff", width = 12, height = 6, unit= "in", res = 300)
p.top20.corecrt
dev.off()



#################################
######5-Core_CRT contribution#########
################################

###calculation of the contribution of core and crt microbiome (all sites)
#1. data preparation
###CRT##
a.crt <- asv.crt.n[,names(asv.crt.n) %in% rownames(mm.a)] #asv of CRT
a.crt <- a.crt[rowSums(a.crt)>0,] #2677

m.crt <- asv.crt.n[,names(asv.crt.n) %in% rownames(mm.m)] #asv of CRT
m.crt <- m.crt[rowSums(m.crt)>0,] #1896

#agri_ditch_control
ac.crt <- asv.crt.n[,names(asv.crt.n) %in% rownames(mm.ac)] #asv of CRT
ac.crt <- ac.crt[rowSums(ac.crt)>0,] #1837
#agri_ditch_treated
am.crt <- asv.crt.n[,names(asv.crt.n) %in% rownames(mm.am)] #asv of CRT
am.crt <- am.crt[rowSums(am.crt)>0,] #1767



####COre microbiome 
a.core <- asv.core[,names(asv.core) %in% rownames(mm.a)] #asv of CRT
a.core <- a.core[rowSums(a.core)>0,] #2677

m.core <- asv.core[,names(asv.core) %in% rownames(mm.m)] #asv of CRT
m.core <- m.core[rowSums(m.core)>0,] #1896

#agri_ditch_control
ac.core <- asv.core[,names(asv.core) %in% rownames(mm.ac)] #asv of CRT
ac.core <- ac.core[rowSums(ac.core)>0,] #1837
#agri_ditch_treated
am.core <- asv.core[,names(asv.core) %in% rownames(mm.am)] #asv of CRT
am.core <- am.core[rowSums(am.core)>0,] #1767

a <- otu.am 
a.m <- mm.am 
xx.crt <- am.crt

levels(factor(a.m$YW))
counts <- table(factor(a.m$YW)) #count sample numbers in each level

# Step 2: Keep levels with more than one sample
valid_levels <- names(counts[counts > 1])#58 for agr; 71 for mixed;73 for all samples


###Calculation
nReads <- 1018
#library(vegan)

# Initialize output

out_df1 <- data.frame()
for ( date in valid_levels) { 
  names <-rownames(a.m)[as.character(a.m$YW) %in% date]
  a.d <- a[,names]
  a.d <- a.d[rowSums(a.d)>0,]
  a.d.m <- a.m[a.m$YW == date,]
  dim(a.d)
  dim(a.d.m)
  
  
  # calculating BC dissimilarity based on the 1st ranked OTU
  
  x <- apply(combn(ncol(a.d), 2), 2, function(x) sum(abs(a.d[,x[1]]- a.d[,x[2]]))/(2*nReads))
  
  a.crt.d <- xx.crt[,names]
  a.crt.d <- a.crt.d[rowSums(a.crt.d)>0,]
  dim(a.crt.d)
  
  x.crt <- apply(combn(ncol(a.crt.d), 2), 2, function(x) sum(abs(a.crt.d[,x[1]]- a.crt.d[,x[2]]))/(2*nReads))
  out_df1 <- rbind(out_df1, data.frame(
    Date=date,
    dist_all=sum(x),
    dist_core=sum(x.crt),
    contribution = (sum(x.crt)/sum(x)*100)))
}



write.csv(out_df1, "allsites_crt_contribution_agr_20250628.csv", quote=FALSE)
write.csv(out_df1, "allsites_crt_contribution_m_20250628.csv", quote=FALSE)
write.csv(out_df1, "allsites_crt_contribution_agr_ac_20250628.csv", quote=FALSE)
write.csv(out_df1, "allsites_crt_contribution_agr_am_20250628.csv", quote=FALSE)
write.csv(out_df1, "allsites_core_contribution_agr_20250628.csv", quote=FALSE)
write.csv(out_df1, "allsites_core_contribution_m_20250628.csv", quote=FALSE)
write.csv(out_df1, "allsites_core_contribution_agr_ac_20250628.csv", quote=FALSE)
write.csv(out_df1, "allsites_core_contribution_agr_am_20250628.csv", quote=FALSE)

##############################################
###########6-Area plot#########
#################################

#area plot
# Change area plot fill colors by groups
df <- read.csv(file="core_crt_contribution_allsites_20250628.csv",sep = ",")
df$Week <- as.factor(df$Week)
df$Year <- as.factor(df$Year)
df$site_type <- as.factor(df$site_type)
df <- df[,-2]

library(tidyr)

df_long <- df %>%
  pivot_longer(
    cols = c(Core, CRT, non_Core.CRT),     # columns to pivot
    names_to = "Fraction",                 # name of new column for original column names
    values_to = "Abundance"                # name of new column for values
  )

names(df_long)
df_long$Fraction <- as.factor(df_long$Fraction)
df_long$Week <- as.numeric(df_long$Week)
#df$Fraction <- as.factor(df$Fraction)

df.a <- df_long[df_long$site_type == "Agriculture", ]
dim(df.a)

df.m <- df_long[grepl("Mixed", df_long$site_type),]
df.ac <- df_long[grepl("Agriculture_Control", df_long$site_type),]
df.am <- df_long[grepl("Agriculture_Management", df_long$site_type),]

str(df_long)
str(df.a)
dim(df.a)
dim(df.m)



p.a<-ggplot(df.a, aes(x=Week, y= Abundance, fill= Fraction)) +
  geom_area() +
  facet_wrap(Year ~ ., nrow = 1)+
  ylab("Contribution (%)") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        strip.text = element_text(size=10),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_blank())+
  ggtitle("agr")
p.a

p.m<-ggplot(df.m, aes(x=Week, y= Abundance, fill= Fraction)) +
  geom_area() +
  facet_wrap(Year ~ ., nrow = 1)+
  ylab("Contribution (%)") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        strip.text = element_text(size=10),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_blank())+
  ggtitle("mixed")
p.m

p.am <- plot_grid(p.a , 
                  p.m, 
                  nrow = 2)
p.ac<-ggplot(df.ac, aes(x=Week, y= Abundance, fill= Fraction)) +
  geom_area() +
  facet_wrap(Year ~ ., nrow = 1)+
  ylab("Contribution (%)") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        strip.text = element_text(size=10),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_blank())+
  ggtitle("(A) agr_control")
p.ac

p.am<-ggplot(df.am, aes(x=Week, y= Abundance, fill= Fraction)) +
  geom_area() +
  facet_wrap(Year ~ ., nrow = 1)+
  ylab("Contribution (%)") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        strip.text = element_text(size=10),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_blank())+
  ggtitle("(B) agr_managed")
p.am

pdf("CRTCore_contribution_allsites_allsamples_20250628.pdf", width = 8, height = 5)
#library(cowplot)
plot_grid(p.a , 
          p.m, 
          nrow = 2)
#print(p)
dev.off()


pdf("Fig.2D_CRTCore_contribution_allsites_allsamples_20250628_nolegend_new.pdf", width = 8, height = 5)
#library(cowplot)
p.am
#print(p)
dev.off()

tiff("Fig.2D_CRTCore_contribution_allsites_allsamples_20250628_nolegend_new.tiff", width = 8, height = 5, unit = "in", res = 300)
#library(cowplot)
plot_grid(p.a + theme(legend.position = "none"), 
          p.m+ theme(legend.position = "none"), 
          nrow = 2)
#print(p)
dev.off()

pdf("supp.FigS5_CRTCore_contribution_allsites_agr_management_20250628_nolegend.pdf", width = 8, height = 5)
#library(cowplot)
plot_grid(p.ac + theme(legend.position = "none"), 
          p.am+ theme(legend.position = "none"), 
          nrow = 2)
#print(p)
dev.off()

tiff("supp.FigS5_CRTCore_contribution_allsites_agr_management_20250628_nolegend.tiff", width = 8, height = 5, unit = "in", res = 300)
#library(cowplot)
plot_grid(p.ac + theme(legend.position = "none"), 
          p.am+ theme(legend.position = "none"), 
          nrow = 2)
#print(p)
dev.off()


#####ANOVA
df2 <- read.csv(file="core_crt_contribution_allsites_20250628_2.csv",sep = ",")
df2$Year <- as.factor(df2$Year)
df2$Week <- as.factor(df2$Week)
df2$site_type <- as.factor(df2$site_type)
df.am <- df2[df2$site_type %in% c("Agriculture", "Mixed"), ]
df.m <- df.am[df.am$site_type == "Mixed",]
df.a <- df.am[df.am$site_type == "Agriculture",]
df.a.s <- df2[df2$site_type %in% c("Agriculture_Treatment","Agriculture_Control"),]
df.a.m <- df.a.s[df.a.s$site_type == "Agriculture_Treatment",]
df.a.nm <- df.a.s[df.a.s$site_type == "Agriculture_Control",]
dim(df.a.m)
str(df.am)

library(lme4)
library(lmerTest)
 mean(df.a$Core)
 mean(df.m$Core)
 mean(df.am$Core)
library(glmmTMB)
 df.a.s$s_y <- paste0(df.a.s$site_type,"_",df.a.s$Year)
model_gamma <- glm(CRT ~ s_y , data = df.a.s, family = Gamma(link = "log"))
summary(model_gamma)

model_gamma <- glm(CRT ~ Year, data = df.a.m, family = Gamma(link = "log"))
summary(model_gamma)

sort(df.a$Core)
library(emmeans)
emm <- emmeans(model_gamma, pairwise ~ Year, adjust = "tukey")
summary(emm$contrasts)



############
##7-Random forest
#################


library(dplyr)
mm.499.water <- mm.499[,water_properties]
mm.499.water <- cbind(mm.499[,c(8,9,151)], mm.499.water)


mm.499.water <- mm.499.water %>%
  mutate(across(-c(1:3), ~ as.numeric(trimws(.))))

df_summary <- mm.499.water %>%
  group_by(site_type, Week, Year) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

df_summary_a <- df_summary %>%
  filter(site_type == "agr")

df_summary_m <- df_summary %>%
  filter(site_type == "mixed")


library(dplyr)

# Convert Year and Week to character in both data frames (to be consistent)
df_summary_m <- df_summary_m %>%
  mutate(
    Year = as.character(Year),
    Week = as.character(Week)
  )

df.m <- df.m %>%
  mutate(
    Year = as.character(Year),
    Week = as.character(Week)
  )
df_combined <- inner_join(df_summary_m, df.m, by = c("Year", "Week"))
dim(df_combined)
names(df_combined)
df_combined_a <- df_combined[,-c(1,2,3,20,21,22)]
df_combined_m <- df_combined[,-c(1,2,3,20,21,22)]

library(randomForest)
df_a_core <- df_combined_a[,1:17]
df_a_crt <- df_combined_a[,c(1:16,18)]
df_m_core <- df_combined_m[,1:17]
df_m_core <- na.omit(df_m_core)
df_m_crt <- df_combined_m[,c(1:16,18)]
df_m_crt <- na.omit(df_m_crt)
set.seed(123)
rf_model.core.a <- randomForest(Core ~ ., data = df_a_core, importance = TRUE, ntree = 500)
rf_model.crt.a <- randomForest(CRT ~ ., data = df_a_crt, importance = TRUE, ntree = 500)
rf_model.core.m <- randomForest(Core ~ ., data = df_m_core, importance = TRUE, ntree = 500)
rf_model.crt.m <- randomForest(CRT ~ ., data = df_m_crt, importance = TRUE, ntree = 500)
print(rf_model.core.a) # R2 = 57.64%
print(rf_model.crt.a)  # R2 = 37.31%
print(rf_model.core.m) # R2 = 22.6%
print(rf_model.crt.m)  # R2 = 7.88%




agr_core_imp = rf_model.core.a$importance %>%
  as.data.frame() %>%
  dplyr::mutate(site_type = 'agr') %>%
  mutate(Fraction = "Core")
agr_core_imp$variables = rownames(agr_core_imp)

agr_crt_imp = rf_model.crt.a$importance %>%
  as.data.frame() %>%
  dplyr::mutate(site_type = 'agr') %>%
  mutate(Fraction = "CRT")
agr_crt_imp$variables = rownames(agr_crt_imp)

m_core_imp = rf_model.core.m$importance %>%
  as.data.frame() %>%
  dplyr::mutate(site_type = 'mixed') %>%
  mutate(Fraction = "Core")
m_core_imp$variables = rownames(m_core_imp)

m_crt_imp = rf_model.crt.m$importance %>%
  as.data.frame() %>%
  dplyr::mutate(site_type = 'mixed') %>%
  mutate(Fraction = "CRT")
m_crt_imp$variables = rownames(m_crt_imp)



# variable naming order to match that from agriculture (set factor levels)

agr_core_imp = agr_core_imp %>% 
  dplyr::arrange((`%IncMSE`))
agr_core_imp$variables = factor(agr_core_imp$variables, levels = agr_core_imp$variables)
agr_crt_imp$variables = factor(agr_crt_imp$variables, levels = agr_core_imp$variables)
m_core_imp$variables = factor(m_core_imp$variables, levels = agr_core_imp$variables)
m_crt_imp$variables = factor(m_crt_imp$variables, levels = agr_core_imp$variables)


df_plot = rbind(agr_core_imp, agr_crt_imp, m_core_imp,m_crt_imp) %>%
  dplyr::rename(IncMSE = `%IncMSE`)

df_plot$site_type = factor(df_plot$site_type, levels = c('agr', 'mixed'))
df_plot$Fraction = factor(df_plot$Fraction, levels = c('Core', 'CRT'))



## plot the importance ##
# facet labels

flabs = c('agr', 'mixed')
names(flabs) = c('agr', 'mixed', 'forest')

r2_df <- data.frame(
  site_type = c("agr", "agr", "mixed", "mixed"),
  Fraction = c("Core", "CRT", "Core", "CRT"),
  R2 = c(57.64, 37.31, 22.6, 7.88),  # numeric values only
  variables = "",                   # dummy x for geom_text
  IncMSE = 0                        # dummy y for geom_text
)


p.core.crt.rf = ggplot(data = df_plot, aes(x = variables, y = IncMSE)) +
  geom_segment(mapping = aes(x=variables, xend=variables, y=0, yend=IncMSE), colour = "#0075DC") +
  geom_point(mapping = aes(size = IncNodePurity), colour = "#0075DC", alpha = 0.5) +
  geom_text(
    data = r2_df,
    aes(label = paste0("R² = ", R2, "%")),
    x = -Inf, y = -Inf,
    hjust = -0.5, vjust = -1.2,
    inherit.aes = FALSE
  ) +
  theme_bw() +
  coord_flip() +
  theme(legend.position = 'bottom',
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  facet_grid(site_type ~ Fraction, labeller = labeller(site_type = flabs))+
    theme(
      axis.text = element_text(size = 10, colour = "black"),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold"),
  )




pdf("randomforeat_corecrtcontribution_waterproperties_20250628.pdf", width = 8, height =8)
p.core.crt.rf
dev.off()

tiff("randomforeat_corecrtcontribution_waterproperties_20250628.tiff", width = 8, height =8, unit = "in", res = 300)
p.core.crt.rf
dev.off()


###################8-relative abundance#########


asv.16S.all.fecal
core.list
asv.all.core <- asv.16S.all.fecal[rownames(asv.16S.all.fecal) %in% core.list,]
dim(asv.all.core)
asv.all.crt <- asv.16S.all.fecal[rownames(asv.16S.all.fecal) %in% rownames(asv.crt.n),]
dim(asv.all.crt)

rel.core <- mean(colSums(asv.all.core)/colSums(asv.16S.all.fecal))
rel.crt <- mean(colSums(asv.all.crt)/colSums(asv.16S.all.fecal))


###taxonomy######
tax.16S.all.fecal <- as.data.frame(tax_table(physeq.16S.all.fecal))
dim(tax.16S.all.fecal)
names(tax.16S.all.fecal)

tax.core <- tax.16S.all.fecal[rownames(tax.16S.all.fecal) %in% rownames(asv.all.core),]
identical(rownames(asv.all.core), rownames(tax.core))
asv.tax.core <- cbind(asv.all.core, tax.core)


tax.crt <- tax.16S.all.fecal[rownames(tax.16S.all.fecal) %in% rownames(asv.all.crt),]
identical(rownames(asv.all.crt), rownames(tax.crt))
asv.tax.crt <- cbind(asv.all.crt, tax.crt)


# Step 1: Split abundance and taxonomy
asv.tax <- asv.tax.crt
abundance <- asv.tax[ , 1:(ncol(asv.tax)-7)]  # assuming last 7 cols are taxonomy
taxonomy <- asv.tax[ , (ncol(asv.tax)-6):ncol(asv.tax)]  # taxonomy cols

# Step 2: Add Genus info to abundance matrix
abundance$Genus <- taxonomy$Genus

# Step 3: Sum abundance by Genus across all ASVs
genus_table <- abundance %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  column_to_rownames("Genus") %>%
  t() %>%
  as.data.frame()

# Result: samples as rows, genera as columns
head(genus_table)
dim(genus_table)
genus.core.ra <- decostand(genus_table, "total", MARGIN = 1)
genus.crt.ra <- decostand(genus_table, "total", MARGIN = 1)
dim(genus.crt.ra) #499 275
dim(genus.core.ra)
sort(colMeans(genus.core.ra),decreasing = TRUE)
sort(colMeans(genus.crt.ra),decreasing = TRUE)

xx <- names(genus.crt.ra)
yy <- gsub("_.*","",xx)

names(genus.crt.ra) <- yy
# Transpose, aggregate by column names, then transpose back
genus.crt.agg <- t(rowsum(t(genus.crt.ra), group = colnames(genus.crt.ra)))

dim(genus.crt.agg) #499 246

pathogen.list
genus.crt.p <- genus.crt.agg[,colnames(genus.crt.agg) %in% pathogen.list]
dim(genus.crt.p) #499  53

mean(rowSums(genus.crt.p)) #0.3381068
sort(colMeans(genus.crt.p), decreasing = TRUE)
