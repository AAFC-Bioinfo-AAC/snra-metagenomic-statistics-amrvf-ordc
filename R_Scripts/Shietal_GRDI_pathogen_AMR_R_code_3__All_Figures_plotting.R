


# Load required packages
library(phyloseq)
library(microbiome)    # for `transform`
library(ggplot2)
library(dplyr)
library(reshape2)





############################################################
####################Section 3.1###############
#########################################################

physeq.16S.fecal.499
physeq.16S.fecal.270
physeq.shotgun.fecal.270

#------------------------
#Fig. 1# 
#---------------------------

# Extract taxonomy tables
tax16S <- as.data.frame(tax_table(physeq.16S.fecal.270))
taxSG  <- as.data.frame(tax_table(physeq.shotgun.fecal.270))

##----Fig. 1A------##

phyla_16S <- unique(na.omit(tax16S$Phylum))
phyla_SG  <- unique(na.omit(taxSG$Phylum))
normalize_phylum <- function(x) {
  # Harmonize GTDB vs SILVA naming differences
  x <- gsub("Actinobacteriota", "Actinobacteria", x)
  x <- gsub("Bacteroidota", "Bacteroidetes", x)
  x <- gsub("Verrucomicrobiota", "Verrucomicrobia", x)
  x <- gsub("Chlamydiota", "Chlamydiae", x)
  x <- gsub("Firmicutes_A", "Firmicutes", x)
  x <- gsub("Firmicutes_B", "Firmicutes", x)
  x <- gsub("Campylobacterota", "Epsilonproteobacteria", x)
  x <- gsub("Patescibacteria", "Candidate_phyla_radiation", x)
  x <- gsub("Myxococcota", "Deltaproteobacteria", x)
  x <- gsub("Desulfobacterota", "Deltaproteobacteria", x)
  x <- gsub("Spirochaetota", "Spirochaetes", x)
  x <- gsub("Planctomycetota", "Planctomycetes", x)
  x <- gsub("Nitrospirota", "Nitrospirae", x)
  x <- gsub("Elusimicrobiota", "Elusimicrobia", x)
  x <- gsub("Fibrobacterota", "Fibrobacteres", x)
  x <- gsub("Gemmatimonadota", "Gemmatimonadetes", x)
  x <- gsub("Acidobacteriota", "Acidobacteria", x)
  x <- gsub("Chloroflexota", "Chloroflexi", x)
  x
}
phyla_16S_norm <- normalize_phylum(phyla_16S)
phyla_SG_norm  <- normalize_phylum(phyla_SG)


shared_phyla <- intersect(phyla_16S_norm, phyla_SG_norm)
unique_16S   <- setdiff(phyla_16S_norm, phyla_SG_norm)
unique_SG    <- setdiff(phyla_SG_norm, phyla_16S_norm)

cat("Shared:", length(shared_phyla),
    "\nUnique to 16S:", length(unique_16S),
    "\nUnique to Shotgun:", length(unique_SG))

library(ggvenn)

# Create a list of phyla sets
phyla_list <- list(
  Amplicon = phyla_16S_norm,
  Shotgun  = phyla_SG_norm
)

# Plot Venn diagram
p.vg.p <- ggvenn(
  phyla_list,
  fill_color = c("#66c2a5", "#fc8d62"),
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 4
)+ 
  ggtitle("Phylum-level Overlap") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
p.vg.p


##--Fig. 1B------##
genusa_16S <- unique(na.omit(tax16S$Genus))
genus_SG  <- unique(na.omit(taxSG$Genus))


# Create a list of phyla sets
genus_list <- list(
  Amplicon = genusa_16S,
  Shotgun  = genus_SG
)

# Plot Venn diagram
p.vg.g <- ggvenn(
  genus_list,
  fill_color = c("#66c2a5", "#fc8d62"),
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 4
)+
  ggtitle("Genus-level Overlap") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
p.vg.g



#--Fig.1C, D---#


##physeq.16S.fecal.499
# Agglomerate to Genus level
phy.genus <- tax_glom(physeq.16S.fecal.499, taxrank = "Genus")

# Transform to relative abundance
phy.genus.a <- otu_table(phy.genus)
phy.genus.a <- as.data.frame(as.matrix(phy.genus.a))
genus.names <- as.character(tax_table(phy.genus)[, "Genus"])
rownames(phy.genus.a)  <-  genus.names 

# Suppose phy.genus.ra is a data frame or matrix

# Step 1: create a new column with simplified names
phy.genus.a$taxa <- sub("_.*", "", rownames(phy.genus.a))

# Step 2: sum rows by the new taxa column
phy.genus.collapsed <- aggregate(. ~ taxa, data = phy.genus.a, FUN = sum)

# Step 3: set rownames and remove the temporary column
rownames(phy.genus.collapsed) <- phy.genus.collapsed$taxa
phy.genus.collapsed$taxa <- NULL

phy.genus.collapsed
phy.genus.collapsed <- t(phy.genus.collapsed)
phy.genus.collapsed <- as.data.frame(phy.genus.collapsed)

phy.genus.ra <- sweep(phy.genus.collapsed, 1, rowSums(phy.genus.collapsed), FUN = "/")

# Get top 20 genera
# Compute mean abundance for each genus
top20.1 <- data.frame(
  Genus = colnames(phy.genus.ra),
  Abundance = colMeans(phy.genus.ra, na.rm = TRUE)
) %>%
  arrange(desc(Abundance)) %>%
  slice(1:20)

top20.1

##physeq.16S.fecal.270

# Agglomerate to Genus level
phy.genus.2 <- tax_glom(physeq.16S.fecal.270, taxrank = "Genus")

# Transform to relative abundance
phy.genus.16s.270 <- otu_table(phy.genus.2)
phy.genus.16s.270 <- as.data.frame(as.matrix(phy.genus.16s.270))
genus.names <- as.character(tax_table(phy.genus.2)[, "Genus"])
rownames(phy.genus.16s.270 )  <-  genus.names 

# Suppose phy.genus.ra is a data frame or matrix

# Step 1: create a new column with simplified names
phy.genus.16s.270$taxa <- sub("_.*", "", rownames(phy.genus.16s.270))

# Step 2: sum rows by the new taxa column
phy.genus.collapsed <- aggregate(. ~ taxa, data = phy.genus.16s.270, FUN = sum)

# Step 3: set rownames and remove the temporary column
rownames(phy.genus.collapsed) <- phy.genus.collapsed$taxa
phy.genus.collapsed$taxa <- NULL

phy.genus.collapsed
phy.genus.collapsed <- t(phy.genus.collapsed)
phy.genus.collapsed <- as.data.frame(phy.genus.collapsed)

phy.genus.ra.2 <- sweep(phy.genus.collapsed, 1, rowSums(phy.genus.collapsed), FUN = "/")

top20.2 <- data.frame(
  Genus = colnames(phy.genus.ra.2),
  Abundance = colMeans(phy.genus.ra.2, na.rm = TRUE)
) %>%
  arrange(desc(Abundance)) %>%
  slice(1:20)

top20.2


###physeq.shotgun.fecal.270
# Agglomerate to Genus level
phy.genus.3 <- tax_glom(physeq.shotgun.fecal.270, taxrank = "Genus")

# Transform to relative abundance
phy.genus.sg.270 <- otu_table(phy.genus.3)
phy.genus.sg.270 <- as.data.frame(as.matrix(phy.genus.sg.270))
genus.names <- as.character(tax_table(phy.genus.3)[, "Genus"])
rownames(phy.genus.sg.270 )  <-  genus.names 

phy.genus.sg.270 <- as.data.frame(t(phy.genus.sg.270))



phy.genus.ra.3 <- sweep(phy.genus.sg.270, 1, rowSums(phy.genus.sg.270), FUN = "/")



top20.3 <- data.frame(
  Genus = colnames(phy.genus.ra.3),
  Abundance = colMeans(phy.genus.ra.3, na.rm = TRUE)
) %>%
  arrange(desc(Abundance)) %>%
  slice(1:20)

top20.3

top.in <- intersect(top20.2$Genus, top20.3$Genus)

top20 <- unique(c(top20.1$Genus, top20.2$Genus,top20.3$Genus))
top20 <- c(top20, "Other")

additional_colours <- c(
  "#A9A9FF", "#FFD700", "#228B22", "#FF69B4", "#00CED1", "#CD5C5C",
  "#00FA9A", "#DA70D6", "#7FFF00", "#FF4500", "#20B2AA", "#B22222",
  "#ADFF2F", "#008B8B", "#800080", "#FF6347", "#4682B4", "#B8860B"
)

high_contrast_colours_40 <- c(
  "#0075DC", "#993F00", "#2BCE48", "#4C005C", "#F0A3FF", "#FFCC99", "#808080",
  "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", "#003380", "#FFA405", "#FFA8BB",
  "#426600", "#FF0010", "#5EF1F2", "#00998F", "#740AFF", "#990000", "#FFFF00",
  "#8DD3C7", additional_colours
)

# Ensure the number of colors matches the number of genera

length(high_contrast_colours_40)  # Should be at least 39

# Create a named vector
high_contrast_colours_named <- setNames(high_contrast_colours_40[1:length(top20)], top20)

top20.x <- top20.3$Genus ##change to top20.2 for 16S.270 

df.x <- phy.genus.ra.3 ##change to phy.genus.ra.2
# Filter for top 20 and collapse others as "Other"
df.x.x <- df.x[,names(df.x) %in% top20.x]
df.x.x$Other <- 1-rowSums(df.x.x)

df.x.x <- cbind(df.x.x, mm.270$site_type)
names(df.x.x)[22] <- "site_type"

# Step 1: Gather into long format
df_long <- df.x.x %>%
  pivot_longer(cols = -site_type, names_to = "Genus", values_to = "RelativeAbundance")

# Step 2: Optional – calculate mean relative abundance per site type
df_summary <- df_long %>%
  group_by(site_type, Genus) %>%
  summarise(Mean_Abundance = mean(RelativeAbundance, na.rm = TRUE)) %>%
  ungroup()


# Step 1: Calculate total abundance across all site types
genus_order <- df_summary %>%
  group_by(Genus) %>%
  summarise(Total = sum(Mean_Abundance, na.rm = TRUE)) %>%
  arrange(desc(Total)) %>%
  pull(Genus)
genus_order_reordered <- c(setdiff(genus_order, "Other"), "Other")
# Step 2: Reorder the Genus factor in df_summary
df_summary$Genus <- factor(df_summary$Genus, levels = genus_order_reordered)
df_summary$site_type <- factor(df_summary$site_type, levels = c("agr", "mixed","forest"))

df.aaaa <- df_summary[df_summary$Genus == "Other",]
# Step 3: Plot

##16S_499 samples
p1 <- ggplot(df_summary, aes(x = site_type, y = Mean_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  ylab("Relative Abundance (%)") +
  xlab("Site Type") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = high_contrast_colours_named)+
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 10)) 
p1

##16S_270 samples
p2 <- ggplot(df_summary, aes(x = site_type, y = Mean_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  ylab("Relative Abundance (%)") +
  xlab("Site Type") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = high_contrast_colours_named)+
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 10)) 
p2

##shotgun_270 samples
p3 <- ggplot(df_summary, aes(x = site_type, y = Mean_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  ylab("Relative Abundance (%)") +
  xlab("Site Type") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = high_contrast_colours_named)+
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 10)) 
p3




library(cowplot)

# Combine Phylum and Genus Venn diagrams side by side
p.vg <- plot_grid(
  p.vg.p + ggtitle("(A) Phylum-level") + theme(plot.title = element_text(hjust = 0.5)),
  p.vg.g + ggtitle("(B) Genus-level") + theme(plot.title = element_text(hjust = 0.5)),
  ncol = 2,
  labels = NULL
)

# Combine with your bar plots below in a 3-row layout
p.aggregate <- plot_grid(
  p.vg,
  p2 + ggtitle("(C) Amplicon_270 samples") + theme(plot.title = element_text(hjust = 0.5)),
  p3 + ggtitle("(D) Shotgun_270 samples") + theme(plot.title = element_text(hjust = 0.5)),
  nrow = 3,
  rel_heights = c(3, 4, 4) # adjust heights if needed
)

# Draw the combined figure
p.aggregate


pdf("section3.1_Fig.1_270comparison.pdf", width = 8, height = 12)
p.aggregate
dev.off()


tiff("section3.1_Fig.1_270comparison.tiff", width = 8, height = 12, unit = "in", res = 300)
p.aggregate
dev.off()





################################
#####Section 3.2.1###################
#####################################3
###################


#------------------------------------
#-----Fig. 2 & supplementary Fig. S1-----------------------
#--------------------------------------


##Statistical analysis###

physeq.16S.fecal.499 <- readRDS("physeq.16S.fecal.499.rds")

# Create a pattern from fecal$Species
mm <- as.data.frame(sample_data(physeq.16S.fecal.499))
alpha_diversity <- estimate_richness(physeq.16S.fecal.499)
dim(alpha_diversity)
names(alpha_diversity)
alpha_diversity$True.Simpson <- 1/(1-alpha_diversity$Simpson)
alpha_diversity$True.Shannon <- exp(alpha_diversity$Shannon)
rownames(alpha_diversity)
identical(rownames(mm),rownames(alpha_diversity))
mm.div.16S.fecal <- cbind(mm, alpha_diversity)
mm.div.16S.fecal$Block <- mm.div.16S.fecal$site

mm.div.16S.fecal$Block <- gsub("18", "18.19",mm.div.16S.fecal$Block)
mm.div.16S.fecal$Block <- gsub("19", "18.19",mm.div.16S.fecal$Block)
mm.div.16S.fecal$Block <- gsub("20", "20.21",mm.div.16S.fecal$Block)
mm.div.16S.fecal$Block <- gsub("21", "20.21",mm.div.16S.fecal$Block)

mm.div.16S.fecal
##mm.div.pathogen
library(tidyverse)
mm.lme = mm.div.16S.fecal %>%
  dplyr::mutate(YW = paste(Year, Week, sep = '-'))
mm.lme$YW = as.factor(mm.lme$YW)
mm.lme$site_type <- factor(mm.lme$site_type, levels = c("agr", "mixed", "forest"))
mm.div.16S.fecal
mm.lme$Block <- factor(mm.lme$Block)

##Conduct GLME model

library(glmmTMB)

alpha_div_lme <- glmmTMB(True.Shannon ~ site_type + (1 | Block) + (1 | YW),
                         data = mm.lme, 
                         family = Gamma(link = "log"))
summary(alpha_div_lme)
# Fit the null model (without the fixed effect 'site_type')
null_model <- glmmTMB(True.Shannon ~ (1 | Block) + (1 | YW), data = mm.lme, family = Gamma(link = "log"))
# Perform likelihood ratio test
anova(alpha_div_lme, null_model, test = "Chisq")
# anova P = 0.2005 for effect of site type on simpson TD
# anova P value = 0.0294 for effect of site type on shannon TD, not significant
# anova P value = 0.04136 for effect of site type on Chao1, not significant
# Perform pairwise comparisons with different adjustment methods
emmeans::emmeans(alpha_div_lme, pairwise ~ site_type, adjust = "bonferroni")  # Bonferroni


#library(lme4)

#alpha_div_lme <- lme4::glmer(True.Simpson ~ site_type + (1 | Block) + (1 | YW), 
#              data = mm.lme, 
#            family = Gamma(link = "log")) #If your data is continuous but highly skewed (non-normal), you can use a Gamma family, which works well for strictly positive continuous data.
#summary(alpha_div_lme)
#plot(alpha_div_lme)
# Fit the null model (without the fixed effect 'site_type')
#null_model <- lme4::glmer(Chao1 ~ (1 | Block) + (1 | YW), 
#           data = mm.lme, 
#          family = Gamma(link = "log"))
# Perform likelihood ratio test
#anova(alpha_div_lme, null_model, test = "Chisq") 
# anova P = 0.2005 for effect of site type on simpson TD
# anova P value = 0.0294 for effect of site type on shannon TD, not significant
# anova P value = 0.1363 for effect of site type on Chao1, not significant


# Perform pairwise comparisons with different adjustment methods
emmeans::emmeans(alpha_div_lme, pairwise ~ site_type, adjust = "bonferroni")  # Bonferroni


#Simpson
#$contrasts
#contrast       estimate     SE  df z.ratio p.value
#agr - mixed       0.137 0.0895 Inf   1.527  0.3803
#agr - forest      0.685 0.1270 Inf   5.396  <.0001
#mixed - forest    0.548 0.1258 Inf   4.360  <.0001
#Shannon
#> emmeans(alpha_div_lme, pairwise ~ site_type, adjust = "bonferroni")        

#$contrasts
#contrast       estimate    SE  df z.ratio p.value
#agr - mixed       0.273 0.0854 Inf   3.201  0.0041
#agr - forest      0.739 0.1231 Inf   6.002  <.0001
#mixed - forest    0.466 0.1228 Inf   3.793  0.0004

#Results are given on the log (not the response) scale. 
#P value adjustment: holm method for 3 tests 

##Chao1
#$contrasts
#contrast       estimate    SE  df z.ratio p.value
#agr - mixed       0.270 0.0902 Inf   2.988  0.0084
#agr - forest      0.354 0.1254 Inf   2.820  0.0144
#mixed - forest    0.084 0.1234 Inf   0.681  1.0000

#------------Fig. 2A-Simpson-based true diversity index--------------

# use stat_test to get the dataframe structure for ggpubr pairwise comparisons 
stat_test = ggpubr::compare_means(True.Simpson ~ site_type, data = mm.lme, method = 't.test') %>%
  dplyr::mutate(y.position = c(50, 75, 100))

# manual P value formatting - check results (Tukey adjusted)
stat_test$custom_p_val = c('0.086', '< 0.001', '<0.001') 

p1 = ggpubr::ggboxplot(data = mm.lme, x = 'site_type', y = 'True.Simpson', fill = 'site_type',
                       notch = TRUE) +
  scale_x_discrete(name = 'Site Type', labels = c('agr', 'mixed', 'forest')) +
  scale_y_continuous(name = 'Simpson-based True Diversity Index', limits = c(0, 200)) +  # Set limits
  stat_pvalue_manual(stat_test, label = 'P {custom_p_val}', tip.length = 0.02) +
  scale_fill_manual(name = 'Site Type', values = high_contrast_colours, labels = c('agr', 'mixed', 'forest')) +
  theme(legend.position = 'none') +
  annotate(
    "text", 
    x = 2,  # Horizontal position (e.g., center of Mixed group)
    y = 200,  # Vertical position near the top
    label = "P < 0.001", 
    size = 5,  # Adjust text size
    hjust = 0.5  # Center alignment
  )

p1



#--------supplementary Fig.S1-Shannon-based true diversity index----------

# use stat_test to get the dataframe structure for ggpubr pairwise comparisons 
stat_test = ggpubr::compare_means( True.Shannon ~ site_type, data = mm.lme, method = 't.test') %>%
  dplyr::mutate(y.position = c(200, 300, 400))

# manual P value formatting - check results (Tukey adjusted)
stat_test$custom_p_val = c('<0.001', '< 0.001', '0.003') 

p2 = ggpubr::ggboxplot(data = mm.lme, x = 'site_type', y = 'True.Shannon', fill = 'site_type',
                       notch = TRUE) +
  scale_x_discrete(name = 'Site Type', labels = c('agr', 'mixed', 'forest')) +
  scale_y_continuous(name = 'Shannon-based True Diversity Index', limits = c(0, 600)) +  # Set limits
  stat_pvalue_manual(stat_test, label = 'P {custom_p_val}', tip.length = 0.02) +
  scale_fill_manual(name = 'Site Type', values = high_contrast_colours, labels = c('agr', 'mixed', 'forest')) +
  theme(legend.position = 'none') +
  annotate(
    "text", 
    x = 2,  # Horizontal position (e.g., center of Mixed group)
    y = 500,  # Vertical position near the top
    label = "P < 0.001", 
    size = 5,  # Adjust text size
    hjust = 0.5  # Center alignment
  )

p2

p2_2 <- p2 + 
  annotate(
    size = 3,  # Adjust text size
    hjust = 0.5  # Center alignment
  )

#------supplementary Fig.S1-Chao1------------

# use stat_test to get the dataframe structure for ggpubr pairwise comparisons 
stat_test = ggpubr::compare_means( Chao1 ~ site_type, data = mm.lme, method = 't.test') %>%
  dplyr::mutate(y.position = c(800, 1000, 1200))

# manual P value formatting - check results (Tukey adjusted)
library(ggpubr)
stat_test$custom_p_val = c('<0.001', '0.07', '= 1.00') 

p3 = ggpubr::ggboxplot(data = mm.lme, x = 'site_type', y = 'Chao1', fill = 'site_type',
                       notch = TRUE) +
  scale_x_discrete(name = 'Site Type', labels = c('agr', 'mixed', 'forest')) +
  scale_y_continuous(name = 'Chao1 Index', limits = c(0, 1500)) +  # Set limits
  stat_pvalue_manual(stat_test, label = 'P {custom_p_val}', tip.length = 0.02) +
  scale_fill_manual(name = 'Site Type', values = high_contrast_colours, labels = c('agr', 'mixed', 'forest')) +
  theme(legend.position = 'none') +
  annotate(
    "text", 
    x = 2,  # Horizontal position (e.g., center of Mixed group)
    y = 1400,  # Vertical position near the top
    label = "P = 0.01", 
    size = 5,  # Adjust text size
    hjust = 0.05  # Center alignment
  )

p3



#-------------Supplementary Fig.S1----------

# Specify the path where you want to save the PDF file
pdf_file <- file.path(wkdir, "Figures/New/16S_allsamples_alphadiversity_fecalcommunity_alpha_indices_sitetype_20250627.pdf")

# Open the PDF device with specified dimensions
pdf(file = pdf_file, width = 6, height = 4)

# Plot your graph
# Make sure you have assigned your plot object to 'plot' variable
# Example: print(plot)
cowplot::plot_grid(p1, p3, ncol = 2)

# Close the PDF device
dev.off()

###Save as tiff
# Specify the path where you want to save the TIFF file
tiff_file <- file.path(wkdir, "Figures/New/16S_allsamples_alphadiversity_fecalcommunity_alpha_indices_sitetype_20250627.tiff")

# Open the TIFF device with specified dimensions and resolution
tiff(file = tiff_file, width = 6, height = 4, units = "in", res = 300)

# Plot your graph using cowplot
cowplot::plot_grid(p1,  p3, ncol = 2)

# Close the TIFF device
dev.off()


#---------- Fig. 2B --------------

p4 = mm.lme %>%
  ggplot() +
  geom_point(aes(x = Week, y = True.Shannon, colour = site_type), alpha = 0.2) +
  stat_smooth(aes(x = Week, y = True.Shannon, colour = site_type), method = 'loess') +
  theme_bw() +
  scale_colour_manual(
    name = 'Site Type',
    values = high_contrast_colours,
    labels = c('agr', 'mixed', 'forest')
  ) +
  scale_x_continuous(name = "Week of Year", breaks = c(16, 20, 24, 28, 32, 36, 40, 44)) +
  scale_y_continuous(name = 'Shannon-based True Diversity Index', breaks = c(0, 100,  200,  300,400)) +
  coord_cartesian(xlim = c(16, 46), ylim = c(-1, 400),expand = 0) +
  #theme(legend.position = 'none') +
  facet_wrap(~Year)
p4


# Specify the path where you want to save the PDF file
pdf_file <- file.path(wkdir, "Figures/16S_allsamples_alphadiversity_fecalcommunity_shannon_sitetype_overtime_withlegend_20250528.pdf")

# Open the PDF device with specified dimensions
pdf(file = pdf_file, width = 8, height = 6)

# Plot your graph
# Make sure you have assigned your plot object to 'plot' variable
# Example: print(plot)
print(p4)

# Close the PDF device
dev.off()

# Specify the path where you want to save the TIFF file
tiff_file <- file.path(wkdir, "Figures/16S_allsamples_alphadiversity_fecalcommunity_shannon_sitetype_overtime_withlegend_20250528.tiff")

# Open the TIFF device with specified dimensions and resolution
tiff(file = tiff_file, width = 8, height = 6, units = "in", res = 300)

# Plot your graph using cowplot
print(p4)

# Close the TIFF device
dev.off()



#----Fig.2C----------

#### Random Forest Analysis of Alpha Diversity ####
water_properties = c(
  'AMIA_AMN',
  'CONDUCTIVITY_MSC',
  'DISS_OXYGEN_MGL',
  'DOC',
  'NITRATE',
  'NITRITE',
  'ORP_MV',
  'PH',
  'TEMP_C',
  'TOC',
  'TOTKN',
  'TOTPHO',
  'TURBIDITY_NTU',
  'RU_DISM3S', # include water discharge data at nearby river because this gives information about snow melt or exceptionally heavy precipitation events
  'rain_mm_7d', # rainfall from last 7 days prior to sampling
  'avg_temp_C_7d' # average 4d air temperature
)

## compute RF for all samples ##
rf_data = mm.lme %>% 
    dplyr::select(all_of(water_properties), True.Shannon)
rf_data <- na.omit(rf_data)

set.seed(123)
all_rf = randomForest::randomForest(
  True.Shannon ~ .,
  data = rf_data,
  importance = T,
  proximity = T,
  ntree = 1000,
  oob.times = 100
)
all_rf
## compute RF for agriculture ##
rf_data_agr = mm.lme %>% 
  dplyr::filter(site_type == 'agr') %>%
  dplyr::select(all_of(water_properties), True.Shannon)
rf_data_agr <- na.omit(rf_data_agr)

set.seed(123)
agr_rf = randomForest::randomForest(
  True.Shannon ~ .,
  data = rf_data_agr,
  importance = T,
  proximity = T,
  ntree = 1000,
  oob.times = 100
)

## compute RF for mixed ##
rf_data_mixed = mm.lme %>% 
  dplyr::filter(site_type == 'mixed') %>%
  dplyr::select(all_of(water_properties), True.Shannon)
rf_data_mixed <- na.omit(rf_data_mixed)
set.seed(123)
mixed_rf = randomForest::randomForest(
  True.Shannon ~ .,
  data = rf_data_mixed,
  importance = T,
  proximity = T,
  ntree = 1000,
  oob.times = 100
)

## compute RF for forest ##
rf_data_forest = mm.lme %>% 
  dplyr::filter(site_type == 'forest') %>%
  dplyr::select(all_of(water_properties), True.Shannon)
rf_data_forest <- na.omit(rf_data_forest)
set.seed(123)
forest_rf = randomForest::randomForest(
  True.Shannon ~ .,
  data = rf_data_forest,
  importance = T,
  proximity = T,
  ntree = 1000,
  oob.times = 100
)

# random forest model summaries - all RF computed on 1000 trees, oob.times = 100; all other parameters
agr_rf
#Call:
# randomForest(formula = True.Shannon ~ ., data = rf_data_agr,      importance = T, proximity = T, ntree = 1000, oob.times = 100) 
#Type of random forest: regression
#Number of trees: 1000
#No. of variables tried at each split: 5

#Mean of squared residuals: 2175.957
#% Var explained: 59.1
mixed_rf
#Call:
#  randomForest(formula = True.Shannon ~ ., data = rf_data_mixed,      importance = T, proximity = T, ntree = 1000, oob.times = 100) 
#Type of random forest: regression
#Number of trees: 1000
#No. of variables tried at each split: 5

#Mean of squared residuals: 1617.582
#% Var explained: 41.16

forest_rf

#Call:
#  randomForest(formula = True.Shannon ~ ., data = rf_data_forest,      importance = T, proximity = T, ntree = 1000, oob.times = 100) 
#Type of random forest: regression
#Number of trees: 1000
#No. of variables tried at each split: 5

#Mean of squared residuals: 211.6015
#% Var explained: 7.32

## format and merge data ##
all_imp = all_rf$importance %>%
  as.data.frame() 
all_imp$variables = rownames(all_imp)

agr_imp = agr_rf$importance %>%
  as.data.frame() %>%
  dplyr::mutate(site_type = 'agr')
agr_imp$variables = rownames(agr_imp)

mixed_imp = mixed_rf$importance %>%
  as.data.frame() %>%
  dplyr::mutate(site_type = 'mixed')
mixed_imp$variables = rownames(mixed_imp)

forest_imp = forest_rf$importance %>%
  as.data.frame() %>%
  dplyr::mutate(site_type = 'forest')
forest_imp$variables = rownames(forest_imp)

# variable naming order to match that from agriculture (set factor levels)
all_imp = all_imp %>% 
  dplyr::arrange((`%IncMSE`))
all_imp$variables = factor(all_imp$variables, levels = all_imp$variables)
names(all_imp)[1] <- "IncMSE"


agr_imp = agr_imp %>% 
  dplyr::arrange((`%IncMSE`))
agr_imp$variables = factor(agr_imp$variables, levels = agr_imp$variables)

mixed_imp$variables = factor(mixed_imp$variables, levels = agr_imp$variables)
forest_imp$variables = factor(forest_imp$variables, levels = agr_imp$variables)

df_plot = rbind(agr_imp, mixed_imp) %>%
  dplyr::rename(IncMSE = `%IncMSE`)
df_plot$site_type = factor(df_plot$site_type, levels = c('agr', 'mixed'))

## plot the importance ##
# facet labels

flabs = c('agr', 'mixed')
names(flabs) = c('agr', 'mixed')
p <- ggplot(data = all_imp, aes(x = variables, y = IncMSE)) +
  geom_segment(mapping = aes(x=variables, xend=variables, y=0, yend=IncMSE), colour = "#0075DC") +
  geom_point(mapping = aes(size = IncNodePurity), colour = "#0075DC", alpha = 0.5) +
  theme_bw() +
  coord_flip() +
  theme(legend.position = 'bottom',
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
        strip.text = element_text(size = 14))    

p5 = ggplot(data = df_plot, aes(x = variables, y = IncMSE)) +
  geom_segment(mapping = aes(x=variables, xend=variables, y=0, yend=IncMSE), colour = "#0075DC") +
  geom_point(mapping = aes(size = IncNodePurity), colour = "#0075DC", alpha = 0.5) +
  theme_bw() +
  coord_flip() +
  theme(legend.position = 'bottom',
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
        strip.text = element_text(size = 14))    + #Adjust size as needed)
  facet_grid(~site_type,
             labeller = labeller(site_type = flabs))

# extract tree R^2 value (% variance explained)
annotation = data.frame(rsq = c(agr_rf$rsq[agr_rf$ntree], mixed_rf$rsq[mixed_rf$ntree]),
                        site_type = factor(c('agr', 'mixed'), levels = c('agr', 'mixed')))
annotation$rsq = round(annotation$rsq, 4)

# Add R^2 annotation with correct positioning in each facet
p5 = p5 +
  geom_text(
    data = annotation, 
    aes(x = 5, y = 150, label = paste0("R^2 == ", rsq)),  # Adjust 'x' and 'y' positions based on your plot's range
    parse = TRUE,
    inherit.aes = FALSE,
    hjust = 0,   # Align text to the left
    size = 5     # Increase size for visibility
      )

# Print the plot
p5

pdf("fecal_shannon_rf_20251001.pdf", width = 8, height = 4)
p5
dev.off()

#---------Combining Fig. 2A-C--------------

# Standardize font sizes
p2 <- p2 + theme(text = element_text(size = 10))  # Adjust size as needed
p4 <- p4 + theme(text = element_text(size = 10))
p5 <- p5 + theme(text = element_text(size = 10))
# Combine p2 and p4 into Row 1, with width ratio 1:2
row1 <- plot_grid(p2, p4, labels = c("(A)", "(B)"), rel_widths = c(1, 2))

# p5 as Row 2 with label
row2 <- plot_grid(p5, labels = "(C)")

# Combine rows with height ratio 2:3
final_plot <- plot_grid(row1, row2, ncol = 1, rel_heights = c(2, 3))


# Save the combined plot to a PDF
pdf("Fig.2_16S_allsamples_shannondiversity.pdf", width = 8, height = 8)  # adjust width/height as needed
print(final_plot)
dev.off()

tiff("Fig.2_16S_allsamples_shannondiversity.tiff", width = 7, height = 9, unit = "in", res = 300)  # adjust width/height as needed
print(final_plot)
dev.off()




#-----------------------------------------
#-----------Supplementary Fig. S2---------------
#----------------------------------------

# Spearman is more robust to non-normality
cor.test(rf_data_agr$TEMP_C, rf_data_agr$True.Shannon ,  method = "spearman")
library(pdp)

partial(agr_rf, pred.var = "TEMP_C") %>%
  autoplot()

####Spearman correlation
df <- rf_data_agr
# Select only numeric columns
df_numeric <- df[, sapply(df, is.numeric)]
library(Hmisc)
# Compute correlation and p-values
cor_res <- rcorr(as.matrix(df_numeric), type = "spearman")
cor_matrix <- cor_res$r
p_matrix <- cor_res$P

library(ggplot2)
library(ggpubr)

# Assuming df is your data frame
# Plot True.Shannon vs Air Temperature
p1 <- ggscatter(df, x = "TEMP_C", y = "True.Shannon",
                add = "reg.line", conf.int = TRUE,
                cor.coef = TRUE, cor.method = "spearman",
                cor.coef.size = 5,
                xlab = "Air Temperature (°C)", ylab = "Shannon-based True Diversity Index")

# Plot True.Shannon vs TOC
p2 <- ggscatter(df, x = "TOC", y = "True.Shannon",
                add = "reg.line", conf.int = TRUE,
                cor.coef = TRUE, cor.method = "spearman",
                cor.coef.size = 5,
                xlab = "Total Organic Carbon (mg/L)", ylab = "Shannon-based True Diversity Index")

# Plot True.Shannon vs Nitrate
p3 <- ggscatter(df, x = "NITRATE", y = "True.Shannon",
                add = "reg.line", conf.int = TRUE,
                cor.coef = TRUE, cor.method = "spearman",
                cor.coef.size = 5,
                xlab = "Nitrate (mg/L)", ylab = "Shannon-based True Diversity Index")

# Arrange in one figure
pdf("supp_FigS2_True.Shannon_SpearmanCorrelation.pdf", width = 5, height = 10)
ggarrange(p1, p2, p3,
          nrow = 3,
          labels = c("(A)", "(B)", "(C)"),
          label.x = 0.90,     # Move label horizontally (0 = left, 1 = right)
          label.y = 0.98)     # Move label vertically (0 = bottom, 1 = top)

dev.off()

tiff("supp_FigS2_True.Shannon_SpearmanCorrelation.tiff", width = 5, height = 10, unit = "in", res = 300)
ggarrange(p1, p2, p3,
          nrow = 3,
          labels = c("(A)", "(B)", "(C)"),
          label.x = 0.90,     # Move label horizontally (0 = left, 1 = right)
          label.y = 0.98)     # Move label vertically (0 = bottom, 1 = top)

dev.off()





#----------------------------------------------------------------------------------------
#-----------Supplementary Fig. S3 Barplot-top 20 genera of core and CRT---------------
#-----------------------------------------------------------------------------------------


#------p.top20.core---------

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



#------p.top20.crt---------



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


#------combine two plots---------
p.top20.corecrt <- plot_grid(p.top20.core + ggtitle("(A) Core"),
                             p.top20.crt + ggtitle("(B) CRT"), ncol = 2)

pdf("Supplementary Fig. S3- Top20genera_coreCRT.pdf", width = 12, height = 6)
p.top20.corecrt
dev.off()

tiff("Supplementary Fig. S3- Top20genera_coreCRT.tiff", width = 12, height = 6, unit= "in", res = 300)
p.top20.corecrt
dev.off()










#----------------------------------------------------------------------------------------
#-----------Supplementary Fig. S4 Random forest for core and crt contribution---------------
#-----------------------------------------------------------------------------------------


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



#-----------------------------------------------------------------
#----Fig. 2D core/CRT area plot & supplementary Fig. S5---------------
#---------------------------------------------------

#area plot
# Change area plot fill colors by groups
df <- read.csv(file="corecrtcontribution_allsites_allsamples.csv",sep = ",")
df$Week <- as.factor(df$Week)
df$Year <- as.factor(df$Year)
df$site_type <- as.factor(df$site_type)

#df$Fraction <- as.factor(df$Fraction)

df.a <- df[grepl("Agriculture", df$site_type),]
df.m <- df[grepl("Mixed", df$site_type),]
df.ac <- df[grepl("Agriculture_Control", df$site_type),]
df.am <- df[grepl("Agriculture_Treatment", df$site_type),]
str(df)
dim(df.all)
dim(df.a)
dim(df.m)
#df$date <- as.Date(df$date)
#df$Week <- strftime(df$date, format = "%V")
library(tidyr)
library(dplyr)

df_long.a <- df.a %>%
  pivot_longer(
    cols = c("Core", "CRT", "non.Core.CRT"),
    names_to = "Fraction",
    values_to = "Contribution"
  )
df_long.a$Fraction <- as.factor(df_long.a$Fraction)
df_long.a$Week <- as.numeric(as.character(df_long.a$Week))

df_long.m <- df.m %>%
  pivot_longer(
    cols = c("Core", "CRT", "non.Core.CRT"),
    names_to = "Fraction",
    values_to = "Contribution"
  )
df_long.m$Fraction <- as.factor(df_long.m$Fraction)
df_long.m$Week <- as.numeric(as.character(df_long.m$Week))

df_long.ac <- df.ac %>%
  pivot_longer(
    cols = c("Core", "CRT", "non.Core.CRT"),
    names_to = "Fraction",
    values_to = "Contribution"
  )
df_long.ac$Fraction <- as.factor(df_long.ac$Fraction)
df_long.ac$Week <- as.numeric(as.character(df_long.ac$Week))

df_long.am <- df.am %>%
  pivot_longer(
    cols = c("Core", "CRT", "non.Core.CRT"),
    names_to = "Fraction",
    values_to = "Contribution"
  )
df_long.am$Fraction <- as.factor(df_long.am$Fraction)
df_long.am$Week <- as.numeric(as.character(df_long.am$Week))
df_long.a <- df_long.a[df_long.a$Year != "2016",]
df_long.m <- df_long.m[df_long.m$Year != "2016",]
cc.a <- ggplot(df_long.a, aes(x = Week, y = Contribution, fill = Fraction)) +
  geom_area() +
  facet_wrap(~ Year, nrow = 1) +  # Updated line
  ylab("Contribution (%)") +
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  ggtitle("(A) agr")

cc.m <- ggplot(df_long.m, aes(x = Week, y = Contribution, fill = Fraction)) +
  geom_area() +
  facet_wrap(~ Year, nrow = 1) +  # Updated line
  ylab("Contribution (%)") +
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  ggtitle("(B) mixed")

cc.ac <- ggplot(df_long.ac, aes(x = Week, y = Contribution, fill = Fraction)) +
  geom_area() +
  facet_wrap(~ Year, nrow = 1) +  # Updated line
  ylab("Contribution (%)") +
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  ggtitle("(C) agr_control")

cc.am <- ggplot(df_long.am, aes(x = Week, y = Contribution, fill = Fraction)) +
  geom_area() +
  facet_wrap(~ Year, nrow = 1) +  # Updated line
  ylab("Contribution (%)") +
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  ggtitle("(D) agr_managed")

pdf("areaplot_ad_mix_allsamples_20250627_without legend.pdf", width = 8, height = 10)
plot_grid(cc.a,cc.m,cc.ac,cc.am, nrow = 4)
dev.off()

tiff("areaplot_ad_mix_allsamples_20250627_without legend.tiff", width = 8, height = 10, unit = "in", res = 300)
plot_grid(cc.a,cc.m,cc.ac,cc.am, nrow = 4)
dev.off()



######################################################
#################Section 3.2.2-Fig. 3_pathogen_siteTypeComparison
#########################################################
# --- Amplicon ---
df.16s <- genus.16s.pathogen.270.ra[, shared.pathogen]
rownames(df.16s) <- rownames(genus.16s.pathogen.270.ra)
df.16s$sampleID <- rownames(df.16s)

df_long.16s <- df.16s %>%
  pivot_longer(-sampleID, names_to = "Genus", values_to = "relativeabundance") %>%
  merge(mm.270[, c("sampleID", "site","site_type", "Year", "Month")], by = "sampleID")

df_long.16s.mix <- df_long.16s %>%
  filter(site_type=="mixed")

summary_df_16s <- df_long.16s %>%
  group_by(Genus, site_type, Year, Month) %>%
  summarise(mean_abund = mean(relativeabundance),
            se_abund = sd(relativeabundance) / sqrt(n()),
            .groups = "drop") %>%
  mutate(dataset = "Amplicon")

summary_df_16s.m <- df_long.16s.mix %>%
  group_by(Genus, site, Year, Month) %>%
  summarise(mean_abund = mean(relativeabundance),
            se_abund = sd(relativeabundance) / sqrt(n()),
            .groups = "drop") %>%
  mutate(dataset = "Amplicon")
summary_df_16s.m$YM <- paste0(summary_df_16s.m$Year,"_",summary_df_16s.m$Month)

summary_df_16s.m$YM <- factor(summary_df_16s.m$YM, levels = c("2016_May" ,"2016_Jun", "2016_Jul", "2016_Sep" ,"2016_Oct" ,"2016_Nov", "2017_May" ,"2017_Jun" ,"2017_Jul",
                                                              "2017_Oct", "2017_Nov", "2018_Apr", "2018_May" ,"2018_Jul", "2018_Aug" ,"2018_Oct" ,"2018_Nov", "2019_Apr",
                                                              "2019_May", "2019_Jun", "2019_Jul" ,"2019_Aug", "2019_Sep", "2019_Oct", "2019_Nov", "2020_Oct", "2020_Nov",
                                                              "2021_May", "2021_Jun", "2021_Jul", "2021_Oct"))
summary_df_16s.m$site <- as.factor(summary_df_16s.m$site)
levels(summary_df_16s.m$site)[1] <- "SN_5"
levels(summary_df_16s.m$site)[2] <- "SN_6"
levels(summary_df_16s.m$site)[3] <- "SN_10"
p.m.site <- ggplot(filter(summary_df_16s.m,Genus == "Acinetobacter"),
                   aes(x = YM, y = mean_abund, color = site, group = site)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_abund - se_abund, ymax = mean_abund + se_abund),
                width = 0.2, size = 0.6) +
  #facet_wrap(~ Genus, scales = "free_y", ncol = 1) +
  scale_color_manual(values = high_contrast_colours) +
  labs(x = "Sampling Time",
       y = "Relative Abundance of Acinetobacter",
       color = "Site") +
  # ggtitle("(A) agr") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("Acinetobacter_site_YM.pdf", width = 8, height = 5)
p.m.site
dev.off()

tiff("Acinetobacter_site_YM.tiff", width = 8, height = 5, unit = "in", res = 300)
p.m.site
dev.off()

# --- Shotgun ---
df.sg <- genus.sg.pathogen.ra[, shared.pathogen]
rownames(df.sg) <- rownames(genus.sg.pathogen.ra)
df.sg$sampleID <- rownames(df.sg)

df_long.sg <- df.sg %>%
  pivot_longer(-sampleID, names_to = "Genus", values_to = "relativeabundance") %>%
  merge(mm.270[, c("sampleID", "site_type", "Year", "Month")], by = "sampleID")

summary_df_sg <- df_long.sg %>%
  group_by(Genus, site_type, Year, Month) %>%
  summarise(mean_abund = mean(relativeabundance),
            se_abund = sd(relativeabundance) / sqrt(n()),
            .groups = "drop") %>%
  mutate(dataset = "Shotgun")

summary_df_all <- bind_rows(summary_df_16s, summary_df_sg)

summary_df_all <- summary_df_all %>%
  mutate(
    site_type = factor(site_type, levels = c("agr", "mixed", "forest")),
    dataset = factor(dataset, levels = c("Amplicon", "Shotgun")),
    Month = factor(Month, levels = month_levels, ordered = TRUE),
    YM = paste0(Year, "_", Month)
  ) %>%
  arrange(Year, Month) %>%
  mutate(YM = factor(YM, levels = unique(YM)))
summary_df_all$Genus <- factor(summary_df_all$Genus, 
                               levels = shared.pathogen, 
                               ordered = TRUE)

p.agr <- ggplot(filter(summary_df_all, site_type == "agr"),
                aes(x = YM, y = mean_abund, color = dataset, group = dataset)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_abund - se_abund, ymax = mean_abund + se_abund),
                width = 0.2, size = 0.6) +
  facet_wrap(~ Genus, scales = "free_y", ncol = 1) +
  scale_color_manual(values = high_contrast_colours) +
  labs(x = "Sampling Time",
       y = "Relative Abundance",
       color = "Dataset") +
  ggtitle("(A) agr") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.mixed <- ggplot(filter(summary_df_all, site_type == "mixed"),
                  aes(x = YM, y = mean_abund, color = dataset, group = dataset)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_abund - se_abund, ymax = mean_abund + se_abund),
                width = 0.2, size = 0.6) +
  facet_wrap(~ Genus, scales = "free_y", ncol = 1) +
  scale_color_manual(values = high_contrast_colours) +
  labs(x = "Sampling Time",
       y = "Relative Abundance",
       color = "Dataset") +
  ggtitle("(B) mixed") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p.mixed

p.forest <- ggplot(filter(summary_df_all, site_type == "forest"),
                   aes(x = YM, y = mean_abund, color = dataset, group = dataset)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_abund - se_abund, ymax = mean_abund + se_abund),
                width = 0.2, size = 0.6) +
  facet_wrap(~ Genus, scales = "free_y", ncol = 1) +
  scale_color_manual(values = high_contrast_colours) +
  labs(x = "Sampling Time",
       y = "Relative Abundance",
       color = "Dataset") +
  ggtitle("(C) forest") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p.forest

p.site.type.comparison <- cowplot::plot_grid(p.agr, p.mixed, p.forest, ncol = 3)

pdf("section3.3_Fig.pathogen_siteTypeComparison.pdf", width = 24, height = 12)
p.site.type.comparison
dev.off()

tiff("section3.3_Fig.pathogen_siteTypeComparison.tiff", width = 24, height = 12, units = "in", res = 300)
p.site.type.comparison
dev.off()



#----------------------
###glmmtesting##
#---------------------
###
#16S
physeq.16S.pathogen.270

#shotgun
physeq.shotgun.pathogen.270


# Step 1: Agglomerate to Genus
phy.genus <- tax_glom(physeq.16S.pathogen.270, taxrank = "Genus")

class(phy.genus)
otu <- as.data.frame(as.matrix(t(otu_table(phy.genus))))

# Step 4: Collapse OTUs manually by cleaned genus
tax.df <- as.data.frame(tax_table(phy.genus))
names(otu) <- tax.df$Genus
otu <- otu[,colSums(otu)>0]
otu.clr <- decostand(otu, "rclr")
dim(otu.clr)
names(otu.clr) <- names(otu)

g.16s.clr <- otu.clr

g.16s.clr.10 <- otu.clr[,names(otu.clr) %in% shared.pathogen]

g.sg.clr.10 <- otu.clr[,names(otu.clr) %in% shared.pathogen]
names(g.16s.clr.10) <- paste0(names(g.16s.clr.10), ".amplicon")
names(g.sg.clr.10) <- paste0(names(g.sg.clr.10), ".shotgun")

g.clr.10 <- cbind(g.16s.clr.10,g.sg.clr.10)
rownames(g.clr.10) <- rownames(otu)

g.clr.10.mm  <- cbind(
  g.clr.10,
  mm.270[, c("site_type", "Year", "Month","sampleID")]
)
g.clr.10.mm$YW <- paste0(g.clr.10.mm$Year,"_",g.clr.10.mm$Month)



#########GLM model############

library(lmerTest)

dim(g.clr.10.mm)
shapiro.test(g.clr.10.mm$Clostridium.amplicon)
model <- glmmTMB(Massilia.shotgun ~ site_type + (1 | YW), data = g.clr.10.mm)
summary(model)
null_model <- glmmTMB(Massilia.shotgun ~ 1 + (1 | YW), 
                      data = g.clr.10.mm)

anova(null_model, model)
#Sphingomonas.amplicon P < 0.001
#Massilia.amplicon P < 0.001
#Microbacterium.amplicon P = 0.52
#Pseudomonas.amplicon P = 0.001
#Clostridium.amplicon P < 0.001
#Clostridium.shotgun P < 0.001


# Estimated marginal means
emm <- emmeans(model, ~ site_type)
# Residual plot
plot(residuals(model) ~ fitted(model))
# Shapiro-Wilk test for residuals
shapiro.test(residuals(model))

# Pairwise comparisons with Tukey adjustment
pairs(emm, adjust = "holm")

library(dplyr)
library(lme4)
library(lmerTest)

metadata_cols <- c("site_type", "Year", "Month", "sampleID", "YW")
abundance_cols <- setdiff(colnames(g.clr.10.mm), metadata_cols)

results <- data.frame(genus = abundance_cols, p_value = NA)
pairwise_results <- list()
for (genus in abundance_cols) {
  df <- g.clr.10.mm %>%
    select(all_of(genus), site_type, YW) %>%
    rename(abundance = all_of(genus)) %>%
    mutate(site_type = factor(site_type), YW = factor(YW))
  
  model <- lmer(abundance ~ site_type + (1 | YW), data = df)
  results$p_value[results$genus == genus] <- anova(model)$`Pr(>F)`[1]  # overall site_type effect
  # Store results in a list
  emm <- emmeans(model, ~ site_type)
  pairwise_results[[genus]] <- pairs(emm, adjust = "tukey")
  
}

results
pairwise_results

write.csv(mm.270, "SNR.meta.270.csv")
library(dplyr)

# Assume g.clr.10.mm is a data.frame
# First, separate metadata and abundance columns
metadata_cols <- c("site_type", "Year", "Month", "sampleID", "YW")
abundance_cols <- setdiff(colnames(g.clr.10.mm), metadata_cols)

# Compute mean abundance for each genus grouped by site_type and YW
g.clr.mean <- g.clr.10.mm %>%
  group_by(site_type, YW) %>%
  summarise(across(all_of(abundance_cols), mean, na.rm = TRUE)) %>%
  ungroup()

# Check
head(g.clr.mean)

# Compute CV per genus for each site_type
genera_cv <- g.clr.10.mm %>%
  group_by(site_type) %>%
  summarise(across(
    all_of(abundance_cols),
    ~ sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE),
    .names = "CV_{.col}"
  )) %>%
  ungroup()

# View result
genera_cv

write.csv(genera_cv, "pathogen.genera.top10.comparison.16s.sg.csv")

model <- glmmTMB(Flavobacterium.amplicon ~ Flavobacterium.shotgun + (1|site_type) + (1|YW), data = g.clr.10.mm)
summary(model) # P < 0.001 E = 0.97

model <- glmmTMB(Legionella.amplicon ~ Legionella.shotgun + (1|site_type) + (1|YW), data = g.clr.10.mm)
summary(model) # P < 0.001, E = 1.22

model <- glmmTMB(Sphingomonas.amplicon ~ Sphingomonas.shotgun + (1|site_type)+ (1|YW) , data = g.clr.10.mm)
summary(model) # P < 0.001, E = 0.28

model <- glmmTMB(Massilia.amplicon ~ Massilia.shotgun + (1|site_type) + (1|YW) , data = g.clr.10.mm)
summary(model) # P = 0.003, E = 0.26

model <- glmmTMB(Pseudomonas.amplicon ~ Pseudomonas.shotgun + (1|site_type) + (1|YW), data = g.clr.10.mm)
summary(model) # P < 0.001, E = 0.439

model <- glmmTMB(Clostridium.amplicon ~ Clostridium.shotgun + (1|site_type) + (1|YW), data = g.clr.10.mm)
summary(model) # P = 0.468, E < 0.01

model <- glmmTMB(Microbacterium.amplicon ~ Microbacterium.shotgun + (1|site_type) + (1|YW), data = g.clr.10.mm)
summary(model) # P = 0.008, E = -0.08









################################################
##########Section 3.3-Fig. 4##################
##########################################

#######AMR

library(tidyr)
library(dplyr)
amr <- read.csv("MG_AMR_read_mapping_80_cutoff_nosnp 2.csv", header = TRUE, sep = ",")

##order the amr list according to the importance
amr.list <- c("imp","oxa","sme",
              "ant3-Dprime", "ant6","ant9",
              "aph3-Dprime",
              "aph6",
              "dfrB",
              "ermF",
              "lnuC",
              "lsa",
              "lsae",
              "mphG",
              "mefA",
              "mefB",
              "mefC",
              "cat",
              "tet39",
              "tet44", 
              "tetM",
              "tetQ",
              "tetT", 
              "tetW",
              "emrAsm",
              "emrBsm",
              "mexF",
              "sodA",
              "sodB",
              "merC",
              "merF",
              "merP",
              "merT",
              "ruvB"
)




#--------------------------------
#----Fig.4A,B------------
#----------------------------





##Convert to presence table
amr1 <- amr %>%
  group_by(sampleID) %>%
  summarize(Count = sum(presence), .groups = "drop")

amr <- amr %>%
  mutate(Site_Year = paste(Site, Year, sep = "_"))

##summary table
summary_pa <- amr %>%
  group_by(AMR_gene, Site_Year) %>%
  summarize(Present = max(presence), .groups = "drop")
summary_pa$AMR_gene <- factor(summary_pa$AMR_gene, levels = amr.list)

# Step 1: Create a lookup table of AMR_gene and Category
amr_category_lookup <- amr %>%
  dplyr::select(AMR_gene, Category) %>%
  distinct()


# Step 2: Join the category info to summary_pa
summary_pa <- summary_pa %>%
  left_join(amr_category_lookup, by = "AMR_gene")
summary_pa$Site_Year <- factor(summary_pa$Site_Year, levels = c("SN_5_2016","SN_5_2017","SN_5_2018","SN_5_2019","SN_5_2020",
                                                                "SN_6_2020","SN_6_2021","SN_10_2017","SN_10_2021","SN_18_2021",
                                                                "SN_19_2017","SN_19_2018","SN_19_2019","SN_19_2020",
                                                                "SN_20_2019","SN_20_2020","SN_21_2018","SN_21_2021","SN_24_2018"))

ggplot(summary_pa, aes(x = Site_Year, y = AMR_gene, fill = factor(Present))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "darkblue")) +
  labs(x = "Site-Year", y = "AMR Gene", fill = "Presence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        legend.position = "none")

amr.heatmap <- ggplot(summary_pa, aes(x = Site_Year, y = AMR_gene, fill = factor(Present))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "blue")) +
  labs(x = "Site-Year", y = "AMR Gene", fill = "Presence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        legend.position = "none")


pdf("AMR_site_Year_heatmap.pdf", width = 6, height = 6)
amr.heatmap
dev.off()

tiff("AMR_site_Year_heatmap.tiff", width = 6, height = 6, unit = "in", res = 300)
amr.heatmap
dev.off()


###statistical analysis

library(vegan)

library(dplyr)
library(tidyr)

# Step 1: Ensure your `amr` data has only the required columns
amr_subset <- dplyr::select(amr, sampleID, AMR_gene, presence)



# Step 2: Pivot wider to make AMR_gene as columns
library(reshape2)

amr_wide <- dcast(
  amr_subset,
  sampleID ~ AMR_gene,
  value.var = "presence",
  fill = 0
)


# Step 3: Set rownames as sampleID
amr_matrix <- as.data.frame(amr_wide)
rownames(amr_matrix) <- amr_matrix$sampleID
#amr_matrix$sampleID <- NULL



# Check result
dim(amr_matrix)



# Extract AMR presence/absence matrix

mm.amr <- mm.270[rownames(mm.270) %in% rownames(amr_matrix),]

# Ensure sampleID is a column in mm.270 (not just rownames)
mm.270$sampleID <- rownames(mm.270)

# Merge mm.270 sample list with amr_wide by sampleID (left join keeps all mm.270 samples)
amr_full <- merge(mm.270[, "sampleID", drop = FALSE], amr1, by = "sampleID", all.x = TRUE)

# Replace all NAs (i.e., missing AMR genes) with 0
amr_full[is.na(amr_full)] <- 0
rownames(amr_full) <- amr_full$sampleID
identical(rownames(amr_full), rownames(mm.270))
amr_full <- amr_full[,-1]
amr_full <- as.data.frame(amr_full)
names(amr_full)[1] <- "amr_count"
amr_full$site_type <- mm.270$site_type

library(MASS)
nb_model <- glm.nb(amr_count ~ site_type, data = amr_full)
summary(nb_model)

anova(nb_model, test = "Chisq")

library(emmeans)
emmeans(nb_model, pairwise ~ site_type)
amr_full$site_type <- factor(amr_full$site_type,levels = c("agr","mixed","forest"))


#-----graph--------
amr.p1 <- ggplot(amr_full, aes(x = site_type, y = amr_count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.6) +
  labs(y = "AMR gene presence", x = "Site type") +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

amr.heatmap <- ggplot(summary_pa, aes(x = Site_Year, y = AMR_gene, fill = factor(Present))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "blue")) +
  labs(x = "Site-Year", y = "AMR Gene", fill = "Presence") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        legend.position = "none")
#amr.p1 <- amr.p1 + theme_minimal(base_size = 12)
#amr.heatmap <- amr.heatmap + theme_minimal(base_size = 12)

combined_plot <- plot_grid(
  amr.p1, 
  amr.heatmap, 
  ncol = 2,       # one column = vertical stack
  rel_widths =   c(1, 3),  # gives more height to heatmap
  labels = c("(A)", "(B)")      # auto labels the plots
)

# Display
print(combined_plot)

pdf("Fig3_AB_AMR_gene_sitetype_heatmap.pdf", width = 10, height = 8)
combined_plot
dev.off()

tiff("Fig3_AB_AMR_gene_sitetype_heatmap.tiff", width = 10, height = 8, unit = "in", res = 300)
combined_plot
dev.off()


#--------------------
#--Fig.4C------------
-------------------------
  
  
###cdc/who list
list.i <- c("Pseudomonas",
            "Rickettsia",
            "Acinetobacter",
            "Mycobacterium",
            "Coxiella",
            "Aeromonas",
            "Leptospira",
            "Escherichia",
            "Erysipelothrix",
            "Staphylococcus",
            "Enterococcus",
            "Brucella",
            "Neisseria",
            "Streptococcus",
            "Campylobacter",
            "Borrelia",
            "Bartonella",
            "Yersinia",
            "Salmonella",
            "Capnocytophaga",
            "Haemophilus",
            "Chlamydia",
            "Shigella",
            "Streptobacillus",
            "Arcobacter")
genus.sg.pathogen.ra.top20 <- genus.sg.pathogen.ra[,list.i]
genus.sg.pathogen.ra.mm <- cbind(mm.270[,c(5,8,9)],genus.sg.pathogen.ra.top20)
genus.sg.pathogen.ra.mm$site_Year <- paste0(genus.sg.pathogen.ra.mm$site, "_", genus.sg.pathogen.ra.mm$Year)
#heatmap
heatmap_data <- genus.sg.pathogen.ra.mm %>%
  group_by(site_Year) %>%
  summarise(across(where(is.numeric), mean)) %>%
  as.data.frame()
#write.csv(heatmap_data,"genus_sitetype_year.csv")
# Move rownames and format for heatmap
rownames(heatmap_data) <- heatmap_data$site_Year
names(heatmap_data)
heatmap_data <- heatmap_data[ , -1]  # remove 'Group' column for heatmap
dim(heatmap_data)
# Optional: log transform to reduce dominance of abundant taxa
heatmap_data_log <- log10(heatmap_data + 1e-5)
library(pheatmap)


# Transpose so genera are rows, site_type-Year groups are columns
p.amr.heatmap <-pheatmap(t(heatmap_data),
                         scale = "row",  # z-score within each genus
                         clustering_distance_cols = "euclidean",
                         clustering_distance_rows = "euclidean",
                         clustering_method = "complete",
                         fontsize_row = 10,
                         fontsize_col = 10
)


pdf("heatmap_knownAMRcarrier.pdf", width = 10, height = 10)
p.amr.heatmap
dev.off()

# Suppose heatmap_data is your original matrix/data frame
# melt it to long format
df_long <- melt(as.matrix(heatmap_data))
df_long$RelativeAbundance <- as.numeric(as.character(df_long$RelativeAbundance))
df_long <- df_long[-c(1:46),]
# Rename columns to something meaningful for ggplot
colnames(df_long) <- c("site_Year", "Genus", "RelativeAbundance")
df_long$site_Year <- as.factor(df_long$site_Year)
ggplot(df_long, aes(x = site_Year, y = Genus, fill = RelativeAbundance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +  # or any color scale you prefer
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10)) +
  labs(fill = "Relative Abundance")


pdf("genera contributing to AMR_top.pdf", width = 10, height = 8)
p.amr.heatmap
dev.off()

tiff("genera contributing to AMR_top.tiff", width = 10, height = 8, unit = "in", res = 300)
p.amr.heatmap
dev.off()





##################################
####Supplementary Fig. S6
#######################################

genus.16s.fecal.270 <- genus.16S.fecal[row.names(genus.16S.fecal) %in% row.names(mm.270),]
dim(genus.16s.fecal.270) #[1] 270 947
genus.16s.fecal.270 <- genus.16s.fecal.270[,colSums(genus.16s.fecal.270)>0] #894
genus.16s.fecal.270.ra <- decostand(genus.16s.fecal.270, "total", MARGIN = 1)

names(genus.16s.fecal.270) <- sub("_.*", "", names(genus.16s.fecal.270))
names(genus.16s.fecal.270.ra) <- sub("_.*", "", names(genus.16s.fecal.270.ra))

# Step 1: Identify the unique column names (after potential cleaning)
unique_names <- unique(colnames(genus.16s.fecal.270))

# Step 2: Aggregate columns with the same name
genus.16s.fecal.270_aggregated <- sapply(unique_names, function(name) {
  cols <- which(colnames(genus.16s.fecal.270) == name)
  rowSums(genus.16s.fecal.270[, cols, drop = FALSE])
})

genus.16s.fecal.270.ra_aggregated <- sapply(unique_names, function(name) {
  cols <- which(colnames(genus.16s.fecal.270.ra) == name)
  rowSums(genus.16s.fecal.270.ra[, cols, drop = FALSE])
})



matching_columns_16S <- names(genus.16s.fecal.270)[str_detect(names(genus.16s.fecal.270), str_c(pathogen.list, collapse = "|"))]
genus.16s.pathogen.270 <- genus.16s.fecal.270[,matching_columns_16S]
dim(genus.16s.pathogen.270) #244

genus.16s.pathogen.270.ra <- genus.16s.fecal.270.ra[,matching_columns_16S]
dim(genus.16s.pathogen.270.ra) #244
mean(rowSums(genus.16s.pathogen.270.ra)) #0.2981
genus.16s.pathogen.270.ra.2 <- genus.16s.pathogen.270.ra
names(genus.16s.pathogen.270.ra.2) <- sub("_.*", "", names(genus.16s.pathogen.270.ra.2))
names(genus.16s.pathogen.270) <- sub("_.*", "", names(genus.16s.pathogen.270))
# Step 1: Identify the unique column names (after potential cleaning)
unique_names <- unique(colnames(genus.16s.pathogen.270))

# Step 2: Aggregate columns with the same name
genus.16s.pathogen.270.ra_aggregated <- sapply(unique_names, function(name) {
  cols <- which(colnames(genus.16s.pathogen.270.ra.2) == name)
  rowSums(genus.16s.pathogen.270.ra.2[, cols, drop = FALSE])
})

genus.16s.pathogen.270_aggregated <- sapply(unique_names, function(name) {
  cols <- which(colnames(genus.16s.pathogen.270) == name)
  rowSums(genus.16s.pathogen.270[, cols, drop = FALSE])
}) 
# Step 3: Convert to data frame
genus.16s.pathogen.270.ra_aggregated <- as.data.frame(genus.16s.pathogen.270.ra_aggregated)

dim(genus.16s.pathogen.270.ra_aggregated) ##157

genus.16s.pathogen.270_aggregated <- as.data.frame(genus.16s.pathogen.270_aggregated)


###shotgun sequencing data

genus.fecal <- genus.sg.pathogen
dim(genus.fecal)
genus.fecal <- as.data.frame(as.matrix(genus.fecal))
genus.fecal.ra <- decostand(genus.fecal, "total",MARGIN = 1)
genus.sg.pathogen <-genus.fecal[,names(genus.fecal) %in% pathogen.list]
genus.sg.pathogen.rclr <- decostand(genus.sg.pathogen, "rclr")
dim(genus.sg.pathogen.rclr)
names(genus.sg.pathogen.rclr) <- names(genus.sg.pathogen)
genus.sg.pathogen.rclr <- genus.sg.pathogen.rclr[,colMeans(genus.sg.pathogen.rclr)>0]
dim(genus.sg.pathogen) #270,230
dim(genus.sg.pathogen.rclr)#270 110

genus.sg.pathogen.ra <- genus.fecal.ra[,names(genus.fecal.ra) %in% pathogen.list]
mean(rowSums(genus.sg.pathogen.ra)) #0.3708653

###Venn diagram
install.packages("VennDiagram")

library(VennDiagram)
# Define the two sets
set1 <- names(genus.16s.pathogen.270_aggregated)
set2 <- names(genus.sg.pathogen)

# Create a list
venn_list <- list("16S" = set1, "Shotgun" = set2)

# Plot the Venn diagram
venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,  # Set to NULL to plot in RStudio
  fill = c("skyblue", "orange"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  main = "Venn Diagram of Fecal-Associated Pathogens"
)

grid::grid.draw(venn.plot)

# Install if not already installed
install.packages("ggvenn")

library(ggvenn)

# Prepare the list
venn_data <- list(
  "Metabarcoding" = names(genus.16s.pathogen.270.ra_aggregated),
  "Shotgun" = names(genus.sg.pathogen.ra)
)

# Plot
pdf("venndiagram_16S_shotgun_pathogen.pdf")
ggvenn(venn_data, fill_color = c("skyblue", "orange"), show_percentage = TRUE)
dev.off()

tiff("venndiagram_16S_shotgun_pathogen.tiff")
ggvenn(venn_data, fill_color = c("skyblue", "orange"), show_percentage = TRUE)
dev.off()


###################################################
##############Supplementary Fig. S7_detection frequencing
###################################################

###detection frequency

genus.16s.pathogen.270
genus.sg.pathogen

genus.16s.10 <- genus.16s.pathogen.270[,shared.pathogen]
genus.sg.10 <- genus.sg.pathogen[,shared.pathogen]

head(genus.16s.10)
head(genus.sg.10)

detect_freq_16s <- colSums(genus.16s.10 > 0) / nrow(genus.16s.10) * 100
detect_freq_sg  <- colSums(genus.sg.10  > 0) / nrow(genus.sg.10)  * 100

detect_freq_df <- data.frame(
  Genus = names(detect_freq_16s),
  Amplicon = detect_freq_16s,
  Shotgun = detect_freq_sg
)

detect_freq_df_long <- reshape2::melt(detect_freq_df, id.vars = "Genus", 
                                      variable.name = "Method", value.name = "Detection")
detect_freq_df_long$Genus <- factor(detect_freq_df_long$Genus, levels = shared.pathogen)
p.df <- ggplot(detect_freq_df_long, aes(x = Genus, y = Detection, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(
    y = "Detection frequency (%)", 
    x = NULL
  ) +
  scale_fill_manual(
    values = c("Amplicon" = "blue", "Shotgun" = "red")
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

pdf("section3.3_supp.pathogen.top10.detectionfreq.16s.sg.pdf", width = 6, height = 4)
p.df
dev.off()

tiff("section3.3_supp.pathogen.top10.detectionfreq.16s.sg.tiff", width = 6, height = 4, unit = "in", res = 300)
p.df
dev.off()



###############################################
##########supplementary Fig. S9_DTW_plot
##################################################

df_compare <- summary_df_all %>%
  select(Genus, site_type, YM, dataset, mean_abund) %>%
  pivot_wider(names_from = dataset, values_from = mean_abund) %>%
  drop_na()

library(dtw)
df_compare_norm <- df_compare %>%
  group_by(Genus, site_type) %>%
  mutate(
    Amplicon_z = scale(Amplicon)[,1],
    Shotgun_z = scale(Shotgun)[,1]
  ) %>%
  ungroup()

dtw_results <- df_compare_norm %>%
  group_by(Genus, site_type) %>%
  summarise(
    dtw_distance = tryCatch(
      {
        if (length(Amplicon_z) > 2 && length(Shotgun_z) > 2)
          dtw(Amplicon_z, Shotgun_z, keep.internals = FALSE, step.pattern = symmetric2)$distance
        else
          NA_real_
      },
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

write.csv(dtw_results,"dtw_comparison_amplicon_shotgun.csv")

dtw_results <- dtw_results %>%
  group_by(Genus) %>%
  mutate(dtw_norm = dtw_distance / max(dtw_distance))

ggplot(dtw_results, aes(x = site_type, y = dtw_norm)) +
  geom_boxplot(fill = "lightblue") +
  geom_jitter(width = 0.1, alpha = 0.7) +
  theme_bw() +
  labs(y = "Normalized DTW distance", x = "Site type",
       title = "Amplicon–Shotgun temporal matching by site type")
kruskal.test(dtw_norm ~ site_type, data = dtw_results)

dtw.plot <- ggplot(dtw_results, aes(x = Genus, y = dtw_distance, fill = site_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_bw() +
  labs(y = "DTW distance", x = "Genus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = high_contrast_colours)

pdf("section3.3_supp.fig_DTW_plot.pdf", width = 6, height =4)
dtw.plot
dev.off()

tiff("section3.3_supp.fig_DTW_plot.tiff", width = 6, height =4, unit = "in", res = 300)
dtw.plot
dev.off()





###################################################
##########supplementary Fig. S8, 10, 11###################
##########################################################

act <- read.csv("qPCR_Acinetobacter baumannii.csv", head= TRUE, sep = ",")
lis <- read.csv("qPCR_Listeria.csv", head= TRUE, sep = ",")

act.m270 <- act[act$sampleID %in% mm.270$sampleID,]
lis.m270 <- lis[lis$sampleID %in% mm.270$sampleID,]
dim(lis.m270) #237

###Listeria presence

genus.16s.pathogen.270
genus.sg.pathogen

dim(genus.sg.pathogen)




###Listeria
sg.lis.237 <- genus.sg.pathogen[rownames(genus.sg.pathogen) %in% lis.m270$sampleID,]
ap.lis.237 <- genus.16s.pathogen.270[rownames(genus.16s.pathogen.270) %in% lis.m270$sampleID,]

sg.lis.237$sampleID <- rownames(sg.lis.237)
ap.lis.237$sampleID <- rownames(ap.lis.237)

# Extract only sampleID and Listeria columns from each
sg.df <- sg.lis.237[, c("sampleID", "Listeria")]
ap.df <- ap.lis.237[, c("sampleID", "Listeria")]
pcr.df <- lis.m270[, c("sampleID", "Listeria.Genus","L..monocytogenes")]

sg.s.lis.237 <- s.sg.f.m270[rownames(s.sg.f.m270) %in% sg.lis.237$sampleID,]
sg.s.lis.237$sampleID <- rownames(sg.s.lis.237)
sg.s.df2.lis <- sg.s.lis.237[, c("sampleID", "Listeria_monocytogenes")]
pcr.df2 <- pcr.df[pcr.df$sampleID %in% sg.s.lis.237$sampleID,]
combined.df.s <- sg.s.df2.lis %>%
  left_join(pcr.df2, by = "sampleID")

sg.s.l.df2 <- sg.s.act.252[,c("sampleID", "Listeria")]

# Rename Listeria columns to identify their source
colnames(sg.df)[2] <- "Listeria_sg"
colnames(ap.df)[2] <- "Listeria_ap"
colnames(pcr.df)[2] <- "Listeria_pcr"

# Merge by sampleID
combined.df <- Reduce(function(x, y) merge(x, y, by = "sampleID", all = TRUE),
                      list(sg.df, pcr.df))

# View the result
head(combined.df)

combined.df$Listeria_sg <- ifelse(combined.df$Listeria_sg > 0, 1, 0) #100%
combined.df$Listeria_pcr
percent_presence <- mean(combined.df$Listeria_pcr == 1) *100 #86.9%

combined.df <- combined.df %>%
  left_join(mm.270[, c("sampleID", "Year", "Month")], by = "sampleID")
combined.df.s <- combined.df.s %>%
  left_join(mm.270[, c("sampleID", "Year", "Month")], by = "sampleID")

combined.df$YM <- paste0(combined.df$Month,"_",combined.df$Year)
combined.df$Listeria_sg <- ifelse(combined.df$Listeria_sg  > 0, 1, 0)

names(combined.df)

plot.df <- combined.df %>%
  pivot_longer(cols = c("Listeria_sg", "Listeria_pcr"),
               names_to = "Method",
               values_to = "Presence")

# Make the line plot
ggplot(plot.df, aes(x = YM, y = Presence, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point() +
  theme_minimal() +
  labs(x = "Year-Month", y = "Listeria Presence", color = "Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



###heatmap
combined.df

# Add Listeria_amplicon column (all 0)
combined.df <- combined.df %>%
  mutate(Listeria_amplicon = 0)
names(combined.df)[2] <- "Shotgun"
names(combined.df)[3] <- "qPCR"
names(combined.df)[7] <- "Amplicon"

# Add Listeria_amplicon column (all 0)
combined.df.s <- combined.df.s %>%
  mutate(Listeria_amplicon = 0)
names(combined.df.s)[2] <- "Shotgun"
names(combined.df.s)[4] <- "qPCR"
names(combined.df.s)[7] <- "Amplicon"
combined.df.s$Shotgun <- ifelse(combined.df.s$Shotgun  > 0, 1, 0)
# Pivot to long format for ggplot
long_df <- combined.df %>%
  select(sampleID, Shotgun, qPCR, Amplicon) %>%
  pivot_longer(cols = -sampleID, names_to = "Method", values_to = "Detection") %>%
  mutate(Method = factor(Method, levels = c("Shotgun", "qPCR", "Amplicon")),
         Detection = factor(Detection, levels = c(0,1)))

# Pivot to long format for ggplot
long_df.s <- combined.df.s %>%
  select(sampleID, Shotgun, qPCR, Amplicon) %>%
  pivot_longer(cols = -sampleID, names_to = "Method", values_to = "Detection") %>%
  mutate(Method = factor(Method, levels = c("Shotgun", "qPCR", "Amplicon")),
         Detection = factor(Detection, levels = c(0,1)))

# Plot
heatmap.l <- ggplot(long_df, aes(x = Method, y = sampleID, fill = Detection)) +
  geom_tile(color = "grey70") +
  scale_fill_manual(values = c("0" = "red", "1" = "blue"),
                    labels = c("Not detected", "Detected")) +
  labs(x = "", y = "Sample ID", fill = "", title = "Detection of Listeria") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 5),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
# Plot
heatmap.l.s <- ggplot(long_df.s, aes(x = Method, y = sampleID, fill = Detection)) +
  geom_tile(color = "grey70") +
  scale_fill_manual(values = c("0" = "red", "1" = "blue"),
                    labels = c("Not detected", "Detected")) +
  labs(x = "", y = "Sample ID", fill = "", title = "Detection of Listeria monocytogenes") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 5),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
heatmap.l.s 

percent_presence2 <- mean(combined.df.s$qPCR == 1)


###########Acinetobacter##########


sg.act.252 <- genus.sg.pathogen[rownames(genus.sg.pathogen) %in% act.m270$sampleID,]
ap.act.252 <- genus.16s.pathogen.270[rownames(genus.16s.pathogen.270) %in% act.m270$sampleID,]
sg.act.252$sampleID <- rownames(sg.act.252)
ap.act.252$sampleID <- rownames(ap.act.252)
pcr.act.252 <- act.m270[act.m270$sampleID %in% rownames(sg.act.252),]
# Extract only sampleID and Listeria columns from each
sg.df2 <- sg.act.252[, c("sampleID", "Acinetobacter")]
ap.df2 <- ap.act.252[, c("sampleID", "Acinetobacter")]
pcr.df2 <- pcr.act.252[, c("sampleID", "Acinetobacter.baumannii", "Conc..cells." )]

# Rename Listeria columns to identify their source
colnames(sg.s.df2)[2] <- "Acinetobacter.baumannii_sg"
colnames(ap.df2)[2] <- "Acinetobacter_ap"
colnames(pcr.df2)[2] <- "Acinetobacter.baumannii_pcr"

# Merge by sampleID
combined.df2 <- Reduce(function(x, y) merge(x, y, by = "sampleID", all = TRUE),
                       list(sg.s.df2, ap.df2,pcr.df2))

combined.df2$Acinetobacter.baumannii_sg <- ifelse(combined.df2$Acinetobacter.baumannii_sg > 0, 1, 0) #100%

combined.df2$Acinetobacter_ap <- ifelse(combined.df2$Acinetobacter_ap > 0, 1, 0)


# Collapse OTUs to species level
physeq.Species <- tax_glom(physeq.shotgun.fecal.270, taxrank = "Species")

# Extract OTU table as data frame
s.sg.f.m270 <- as.data.frame(as.matrix(otu_table(physeq.Species)))

# Extract species names as character vector
species.names <- as.character(tax_table(physeq.Species)[, "Species"])



# Assign to rownames
rownames(s.sg.f.m270) <- species.names


s.sg.f.m270 <- t(s.sg.f.m270)
s.sg.f.m270 <- as.data.frame(s.sg.f.m270)


# Collapse OTUs to species level
physeq.Species2 <- tax_glom(physeq.16S.pathogen.270, taxrank = "Species")

# Extract OTU table as data frame
s.16s.p.m270 <- as.data.frame(as.matrix(otu_table(physeq.Species2)))

# Extract species names as character vector
species.names2 <- as.character(tax_table(physeq.Species2)[, "Species"])
# Assign to rownames
rownames(s.16s.p.m270) <- make.unique(species.names2)
s.16s.p.m270 <- t(s.16s.p.m270)
s.16s.p.m270 <- as.data.frame(s.16s.p.m270)


sg.s.act.252 <- s.sg.f.m270[rownames(s.sg.f.m270) %in% sg.act.252$sampleID,]
sg.s.act.252$sampleID <- rownames(sg.s.act.252)
sg.s.df2 <- sg.s.act.252[, c("sampleID", "Acinetobacter_baumannii")]
ap.s.act.252 <- s.16s.p.m270[rownames(s.16s.p.m270) %in% sg.act.252$sampleID,]
ap.s.act.252$sampleID <- rownames(ap.s.act.252)

sg.s.l.df2 <- sg.s.act.252[,c("sampleID", "Listeria")]
ap.df2.2 <- ap.s.act.252[, c("sampleID", "Acinetobacter")]

# Convert to long format
df_long <- combined.df2 %>%
  mutate(Date = as.Date(sub("SN005.", "", sampleID), format = "%Y%m%d")) %>%
  pivot_longer(
    cols = starts_with("Acinetobacter"),
    names_to = "Method",
    values_to = "Presence"
  )

percent_presence2 <- mean(combined.df2$Acinetobacter.baumannii_pcr == 1) *100 #15.89%
percent_presence2 <- mean(combined.df2$Acinetobacter_ap == 1) *100 #28.29457
percent_presence2 <- mean(combined.df2$Acinetobacter.baumannii_sg == 1) *100 #100%

# Pivot to long format for ggplot
long_df2 <- combined.df2 %>%
  pivot_longer(cols = c(Acinetobacter.baumannii_sg, Acinetobacter_ap, Acinetobacter.baumannii_pcr),
               names_to = "Method", values_to = "Detection") %>%
  mutate(Method = factor(Method, levels = c("Acinetobacter.baumannii_sg",
                                            "Acinetobacter.baumannii_pcr",
                                            "Acinetobacter_ap"
  )),
  Detection = factor(Detection, levels = c(0,1)))
levels(long_df2$Method)[1] <- "Shotgun"
levels(long_df2$Method)[2] <- "qPCR"
levels(long_df2$Method)[3] <- "Amplicon"
# Plot
heatmap.ab <- ggplot(long_df2, aes(x = Method, y = sampleID, fill = Detection)) +
  geom_tile(color = "grey70") +
  scale_fill_manual(values = c("0" = "red", "1" = "blue"),
                    labels = c("Not detected", "Detected")) +
  labs(x = "", y = "Sample ID", fill = "", title = "(B) Detection of Acinetobacter baumannii") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 5),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))

library(cowplot)

# Combine the two heatmaps
combined_plot <- plot_grid(
  heatmap.l,
  heatmap.ab,
  ncol = 2,       # 1 column, stacked vertically
  align = "h",    # align vertically
  labels = NULL,  # titles are already in the plots
  rel_heights = c(1, 1)  # adjust relative heights if needed
)

# Display the combined plot
print(combined_plot)

pdf("section3.3_supp._heatmap_three methods comparison of genera detection.pdf", width = 8, height = 16)
combined_plot
dev.off()

tiff("section3.3_supp._heatmap_three methods comparison of genera detection.tiff", width = 8, height = 16, unit = "in", res = 300)
combined_plot
dev.off()
##Fig.S8
# Save to file (optional)
ggsave("Listeria_Acinetobacter_detection_heatmaps.png", combined_plot, width = 8, height = 12, dpi = 300)

pdf("section3.3_supp._heatmap_three methods comparison of Listeria detection.pdf", width = 4, height = 16)
heatmap.l
dev.off()
tiff("section3.3_supp._heatmap_three methods comparison of Listeria detection.tiff", width = 4, height = 16, unit = "in", res = 300)
heatmap.l
dev.off()


########################
##Species identification
######################


combined.df  <- combined.df %>%
  left_join(mm.270[, c("site","sampleID","site_type")], by = "sampleID")


combined.df2 <- combined.df2 %>%
  left_join(mm.270[, c("site","Year","Month","sampleID","site_type")], by = "sampleID")

names(combined.df2)
combined.df2$YM <- paste0(combined.df2$Year,"_",combined.df2$Month)

# Ensure YM is ordered correctly (convert to factor or date)
combined.df2 <- combined.df2 %>%
  mutate(
    YM = factor(YM, levels = unique(YM)),  # preserve chronological order
    site_type = factor(site_type)
  )
combined.df2$site_type <- factor(combined.df2$site_type, levels = c("agr", "mixed", "forest"))
# Plot
# Summarize by site_type and YM
df_summary <- combined.df2 %>%
  group_by(site, YM) %>%
  summarise(
    mean_conc = mean(Conc..cells., na.rm = TRUE),
    se_conc = sd(Conc..cells., na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

##########
#########Fig. S10
####################
# Ensure YM is ordered
df_summary <- df_summary %>%
  mutate(YM = factor(YM, levels = unique(YM)))
df_summary$site_type <- factor(df_summary$site_type, levels = c("agr", "mixed", "forest"))
# Plot
# Plot with SE error bars
p.qpcr <- ggplot(df_summary, aes(x = YM, y = mean_conc, color = site, group = site)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_conc - se_conc, ymax = mean_conc + se_conc), width = 0.2) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = high_contrast_colours) +
  labs(
    x = "Year-Month",
    y = expression(paste("Acinetobacter baumannii (cells/mL) ")),
    color = "Site"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pdf("Acinetobacter.baumannii.qPCR.pdf", height = 5, width = 8)
p.qpcr
dev.off()

#######################
##FIg. S11_plot for A.baumannii..BlaOXA51.gene
########################################
act_plot <- pcr.act.252 %>%
  select(SNR.site, Sampling.date, `A..baumannii..BlaOXA51.gene.`) %>%
  # Optional: convert to presence/absence if desired
  mutate(Presence = ifelse(`A..baumannii..BlaOXA51.gene.` > 0, 1, 0))
act_plot$SNR.site <- as.factor(act_plot$SNR.site)
levels(act_plot$SNR.site) <- paste0("SN_",levels(act_plot$SNR.site))
# Basic heatmap (presence/absence)
p.qpcr.ab.gene <- ggplot(act_plot, aes(x = Sampling.date, y = SNR.site, fill = factor(Presence))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "grey90", "1" = "red"), 
                    labels = c("Not detected", "Detected")) +
  labs(x = "Sampling Date", y = "Site", fill = "A. baumannii\nblaoxa-51") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("qPCR_A.baumannii with blaoxa51gene.heatmap.pdf", width = 8, height = 6)
p.qpcr.ab.gene
dev.off()

tiff("qPCR_A.baumannii with blaoxa51gene.heatmap.tiff", width = 8, height = 6, unit = "in", res = 300)
p.qpcr.ab.gene
dev.off()


