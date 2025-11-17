# Making box plots depicting the difference in physicochemical properties between samples from different land use and season
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1

#set your working directory were everything will be stored and saved
setwd("C:/Users/Sandersonh/OneDrive - AGR-AGR/Desktop/Wen Chen/Project_1_landuse_AMR_virulence/Statistics_Analysis")


#install.packages("tidyverse")
#install.packages(c("dplyr"))

library(dplyr)
library(tidyverse)

metadata <- read.csv("Metadata_clean_AMRVFMG_project.csv", header = TRUE, stringsAsFactors = FALSE)

data_sampleID <- read.csv("table_lme_metadata_AMRVF_taxa_month.csv", header = TRUE, stringsAsFactors = FALSE)

metadata$SAMPLE_ID <- metadata$Sample_ID_Genome

metadata <- metadata %>%
  filter(SAMPLE_ID %in% data_sampleID$SAMPLE_ID)


#install.packages(c("ggplot2", "ggpubr", "gridExtra"))
library(ggplot2)
library(ggpubr)
library(gridExtra)

###Land_use Main Effects###

# Remove missing values
metadata <- metadata[complete.cases(metadata), ]

# Ensure Land_use is a factor
metadata$Land_use <- factor(metadata$Land_use, levels = c("Agri_ditch", "Mixed", "Forested"))
land_use_colors <- c("Agri_ditch" = "black", "Mixed" = "green", "Forested" = "red")
plot_list_LU <- list()

#DOC
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_DOC_LU <- c("0.0024", "0.0056", "0.0018")


# Create boxplot
plot_list_LU[["DOC_LU"]] <- ggplot(metadata, aes(x = Land_use, y = DOC, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$DOC, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_DOC_LU[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$DOC, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_DOC_LU[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$DOC, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_DOC_LU[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = land_use_colors) 

#DOC_MDL
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_DOC_MDL_LU <- c("4.57E-07", "8.21E-05", "9.23E-05")

# Create boxplot
plot_list_LU[["DOC_MDL_LU"]] <- ggplot(metadata, aes(x = Land_use, y = DOC_MDL, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$DOC_MDL, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_DOC_MDL_LU[1])), color = "black",
            size = 6, vjust = -1.0) +
  geom_text(aes(x = 2, y = max(metadata$DOC_MDL, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_DOC_MDL_LU[2])), color = "black",
            size = 6, vjust = -1.0) +
  geom_text(aes(x = 3, y = max(metadata$DOC_MDL, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_DOC_MDL_LU[3])), color = "black",
            size = 6, vjust = -1.0 ) +
  scale_fill_manual(values = land_use_colors)

#PH
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_PH_LU <- c("1.16E-10", "5.32E-07", "2.80E-06")



# Create boxplot
plot_list_LU[["PH_LU"]] <- ggplot(metadata, aes(x = Land_use, y = PH, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$PH, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_PH_LU[1])), color = "black",
            size = 6, vjust = -1.0) +
  geom_text(aes(x = 2, y = max(metadata$PH, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_PH_LU[2])), color = "black",
            size = 6, vjust = -1.0) +
  geom_text(aes(x = 3, y = max(metadata$PH, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_PH_LU[3])), color = "black",
            size = 6, vjust = -1.0 ) +
  scale_fill_manual(values = land_use_colors)

#TOC
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_TOC_LU <- c("0.0025", "0.0055", "0.0017")

# Create boxplot
plot_list_LU[["TOC_LU"]] <- ggplot(metadata, aes(x = Land_use, y = TOC, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$TOC, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOC_LU[1])), color = "black",
            size = 6, vjust = -1.0) +
  geom_text(aes(x = 2, y = max(metadata$TOC, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOC_LU[2])), color = "black",
            size = 6, vjust = -1.0) +
  geom_text(aes(x = 3, y = max(metadata$TOC, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOC_LU[3])), color = "black",
            size = 6, vjust = -1.0 ) +
  scale_fill_manual(values = land_use_colors)

#TOC_MDL
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_TOC_MDL_LU <- c("4.57E-07", "8.21E-05", "9.23E-05")


# Create boxplot
plot_list_LU[["TOC_MDL_LU"]] <- ggplot(metadata, aes(x = Land_use, y = TOC_MDL, fill = Land_use)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$TOC_MDL, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOC_MDL_LU[1])), color = "black",
            size = 6, vjust = -1.0) +
  geom_text(aes(x = 2, y = max(metadata$TOC_MDL, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOC_MDL_LU[2])), color = "black",
            size = 6, vjust = -1.0) +
  geom_text(aes(x = 3, y = max(metadata$TOC_MDL, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOC_MDL_LU[3])), color = "black",
            size = 6, vjust = -1.0 ) +
  scale_fill_manual(values = land_use_colors)


# Arrange the figures in a 2 by 5 grid
grid_row_height <- 20  # Adjust as needed

# Arrange all figures in a 5 by 2 grid with specified height for each row
grid_plot_LU <- arrangeGrob(grobs = plot_list_LU, ncol = 2, nrow = 3)
grid_plot_LU$heights <- unit(rep(grid_row_height, 10), "cm")

# Save the grid as an image (e.g., PNG) with high resolution
ggsave("LME_grid_image_land_use.png", plot = grid_plot_LU, width = 20, height = 100, units = "in", dpi = 300, limitsize = FALSE)

#Season Main Effects
metadata$Season <- factor(metadata$Season, levels = c("Autumn", "Spring", "Summer"))
season_colors <- c("Autumn" = "black", "Spring" = "red", "Summer" = "green")

plot_list_S <- list()
#AMIA_AMN
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_AMIA_AMN_S <- c("0.9766", "1", "0.0496")

# Create boxplot
plot_list_S[["AMIA_AMN_S"]] <- ggplot(metadata, aes(x = Season, y = AMIA_AMN, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$AMIA_AMN, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_AMIA_AMN_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$AMIA_AMN, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_AMIA_AMN_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$AMIA_AMN, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_AMIA_AMN_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)


#avg_temp_C_1d

# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_avgt1d_S <- c("0.0001", "1.3071e-06", "1.9502e-07")

# Create boxplot
plot_list_S[["avgt1d_S"]] <- ggplot(metadata, aes(x = Season, y = avg_temp_C_1d, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$avg_temp_C_1d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_avgt1d_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$avg_temp_C_1d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_avgt1d_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$avg_temp_C_1d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_avgt1d_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)

#avg_temp_C_2d
pairwise_p_values_avgt2d_S <- c("8.6424e-05",
                       "1.3102e-06",
                       "1.5467e-07")
                       
# Create boxplot
plot_list_S[["avgt2d_S"]] <- ggplot(metadata, aes(x = Season, y = avg_temp_C_2d, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$avg_temp_C_2d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_avgt2d_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$avg_temp_C_2d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_avgt2d_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$avg_temp_C_2d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_avgt2d_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)

#avg_temp_C_3d
pairwise_p_values_avgt3d_S <- c("0.0001", "1.2072e-06", "1.7044e-07")

# Create boxplot
plot_list_S[["avgt3d_S"]] <- ggplot(metadata, aes(x = Season, y = avg_temp_C_3d, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$avg_temp_C_3d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_avgt3d_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$avg_temp_C_3d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_avgt3d_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$avg_temp_C_3d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_avgt3d_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)

#avg_temp_C_5d
pairwise_p_values_avgt5d_S <- c("5.4129e-05",
                       "7.6093e-07",
                       "1.2369e-07")
                       
# Create boxplot
plot_list_S[["avgt5d_S"]] <- ggplot(metadata, aes(x = Season, y = avg_temp_C_5d, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$avg_temp_C_5d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_avgt5d_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$avg_temp_C_5d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_avgt5d_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$avg_temp_C_5d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_avgt5d_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)

#avg_temp_C_7d
pairwise_p_values_avgt7d_S <- c("2.2533e-05",
                       "3.1275e-07",
                       "7.1367e-08")
                       
# Create boxplot
plot_list_S[["avgt7d_S"]] <- ggplot(metadata, aes(x = Season, y = avg_temp_C_7d, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$avg_temp_C_7d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_avgt7d_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$avg_temp_C_7d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_avgt7d_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$avg_temp_C_7d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_avgt7d_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)

#BE_DISM3S
pairwise_p_values_BE_DISM3S_S <- c("0.0021",
                       "0.0003",
                       "0.0907")
                       
# Create boxplot
plot_list_S[["BE_DISM3S_S"]] <- ggplot(metadata, aes(x = Season, y = BE_DISM3S, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$BE_DISM3S, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_BE_DISM3S_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$BE_DISM3S, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_BE_DISM3S_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$BE_DISM3S, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_BE_DISM3S_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)

#CONDUCTIVITY_MSC

pairwise_p_values_con_S <- c("0.0037",
                       "0.0080",
                       "0.0041")
# Create boxplot
plot_list_S[["con_S"]] <- ggplot(metadata, aes(x = Season, y = CONDUCTIVITY_MSC, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$CONDUCTIVITY_MSC, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_con_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$CONDUCTIVITY_MSC, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_con_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$CONDUCTIVITY_MSC, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_con_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)

#daily_avg_solar_radiation_wm2
pairwise_p_values_solar_S <- c("0.0007",
                       "6.6951e-06",
                       "2.8321e-06")
                       
# Create boxplot
plot_list_S[["solar_S"]] <- ggplot(metadata, aes(x = Season, y = daily_avg_solar_radiation_wm2, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$daily_avg_solar_radiation_wm2, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_solar_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$daily_avg_solar_radiation_wm2, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_solar_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$daily_avg_solar_radiation_wm2, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_solar_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)

#DISS_OXYGEN_P
pairwise_p_values_DISS_OXY_S <- c("9.8712e-06",
                       "2.8069e-06",
                       "8.6676e-06")
                       

# Create boxplot
plot_list_S[["DISS_OXY_S"]] <- ggplot(metadata, aes(x = Season, y = DISS_OXYGEN_P, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$DISS_OXYGEN_P, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_DISS_OXY_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$DISS_OXYGEN_P, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_DISS_OXY_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$DISS_OXYGEN_P, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_DISS_OXY_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)


#DOC_MDL
pairwise_p_values_DOC_MDL_S <- c("0.0007",
                       "0.0005",
                       "0.0009")
# Create boxplot
plot_list_S[["DOC_MDL_S"]] <- ggplot(metadata, aes(x = Season, y = DOC_MDL, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$DOC_MDL, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_DOC_MDL_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$DOC_MDL, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_DOC_MDL_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$DOC_MDL, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_DOC_MDL_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)

#max_temp_C_1d
pairwise_p_values_maxt_S <- c("1.5267e-05",
                       "6.1844e-07",
                       "1.4364e-07")
                      
# Create boxplot
plot_list_S[["maxt_S"]] <- ggplot(metadata, aes(x = Season, y = max_temp_C_1d, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$max_temp_C_1d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_maxt_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$max_temp_C_1d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_maxt_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$max_temp_C_1d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_maxt_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)

#min_temp_c_1d
pairwise_p_values_mint_S <- c("0.1557",
                       "4.9716e-05",
                       "2.0912e-06")
                       
# Create boxplot
plot_list_S[["mint_S"]] <- ggplot(metadata, aes(x = Season, y = min_temp_C_1d, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$min_temp_C_1d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_mint_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$min_temp_C_1d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_mint_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$min_temp_C_1d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_mint_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)


#NITRATE
pairwise_p_values_nitrate_S <- c("0.0051",
                       "0.0293",
                       "0.1192")

# Create boxplot
plot_list_S[["nitrate_S"]] <- ggplot(metadata, aes(x = Season, y = NITRATE, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$NITRATE, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_nitrate_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$NITRATE, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_nitrate_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$NITRATE, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_nitrate_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)

#PH
pairwise_p_values_PH_S <- c("4.7008e-09",
                       "6.7332e-10",
                       "1.2242e-09")
                      
# Create boxplot
plot_list_S[["PH_S"]] <- ggplot(metadata, aes(x = Season, y = PH, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$PH, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_PH_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$PH, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_PH_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$PH, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_PH_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)


#rain_mm_3d
pairwise_p_values_rain3d_S <- c("0.0017",
                       "0.0024",
                       "0.0004")
                       
# Create boxplot
plot_list_S[["rain3d_S"]] <- ggplot(metadata, aes(x = Season, y = rain_mm_3d, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$rain_mm_3d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_rain3d_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$rain_mm_3d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_rain3d_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$rain_mm_3d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_rain3d_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)

#rain_mm_5d
pairwise_p_values_rain5d_S <- c("0.0015",
                       "0.0007",
                       "0.0009")
                       

# Create boxplot
plot_list_S[["rain5d_S"]] <- ggplot(metadata, aes(x = Season, y = rain_mm_5d, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$rain_mm_5d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_rain5d_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$rain_mm_5d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_rain5d_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$rain_mm_5d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_rain5d_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)

#rain_mm_7d
pairwise_p_values_rain7d_S <- c("0.0009",
                       "0.0005",
                       "0.0006")
                       


# Create boxplot
plot_list_S[["rain7d_S"]] <- ggplot(metadata, aes(x = Season, y = rain_mm_7d, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$rain_mm_7d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_rain7d_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$rain_mm_7d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_rain7d_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$rain_mm_7d, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_rain7d_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)


#TEMP_C
pairwise_p_values_TEMP_S <- c("0.0002",
                       "1.9098e-05",
                       "2.5012e-06")
                       
# Create boxplot
plot_list_S[["TEMP_S"]] <- ggplot(metadata, aes(x = Season, y = TEMP_C, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$TEMP_C, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TEMP_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$TEMP_C, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TEMP_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$TEMP_C, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TEMP_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)


#TOC_F
pairwise_p_values_TOC_S <- c("3.8477e-11",
                       "5.1498e-12",
                       "1.0292e-11")
                       
# Create boxplot
plot_list_S[["TOC_S"]] <- ggplot(metadata, aes(x = Season, y = TOC_F, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$TOC_F, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOC_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$TOC_F, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOC_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$TOC_F, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOC_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)

#TOC_MDL
pairwise_p_values_TOC_MDL_S <- c("0.0007",
                       "0.0005",
                       "0.0009")
                       
# Create boxplot
plot_list_S[["TOC_MDL_S"]] <- ggplot(metadata, aes(x = Season, y = TOC_MDL, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$TOC_MDL, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOC_MDL_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$TOC_MDL, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOC_MDL_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$TOC_MDL, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOC_MDL_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)

#TOTKN
pairwise_p_values_TOTKN_S <- c("0.0113",
                       "0.0205",
                       "0.0012")
                      
# Create boxplot
plot_list_S[["TOTKN_S"]] <- ggplot(metadata, aes(x = Season, y = TOTKN, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$TOTKN, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOTKN_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$TOTKN, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOTKN_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$TOTKN, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOTKN_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)

#TOTPHO
pairwise_p_values_TOTPHO_S <- c("0.2213",
                       "0.3791",
                       "0.0068")
                       

# Create boxplot
plot_list_S[["TOTPHO_S"]] <- ggplot(metadata, aes(x = Season, y = TOTPHO, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$TOTPHO, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOTPHO_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$TOTPHO, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOTPHO_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$TOTPHO, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOTPHO_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)

#TOTPHO_MDL

pairwise_p_values_TOTPHO_MDL_S <- c("0.0129",
                       "0.0172",
                       "0.0012")
                       
# Create boxplot
plot_list_S[["TOTPHO_MDL_S"]] <- ggplot(metadata, aes(x = Season, y = TOTPHO_MDL, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$TOTPHO_MDL, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOTPHO_MDL_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$TOTPHO_MDL, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOTPHO_MDL_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$TOTPHO_MDL, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_TOTPHO_MDL_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)


#TURBIDITY_NTU
pairwise_p_values_turb_S <- c("0.3309",
                       "1",
                       "0.0209")
                       
# Create boxplot
plot_list_S[["turb_S"]] <- ggplot(metadata, aes(x = Season, y = TURBIDITY_NTU, fill = Season)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$TURBIDITY_NTU, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_turb_S[1])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$TURBIDITY_NTU, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_turb_S[2])), color = "black",
            size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$TURBIDITY_NTU, na.rm = TRUE) + 1,
                label = paste("p =", pairwise_p_values_turb_S[3])), color = "black",
            size = 6, vjust = -0.5 ) +
  scale_fill_manual(values = season_colors)

# Arrange the figures in a 2 by 5 grid
grid_row_height <- 20  # Adjust as needed

# Arrange all figures in a 5 by 2 grid with specified height for each row
grid_plot_S <- do.call(grid.arrange, c(grobs = plot_list_S, ncol = 3, nrow = 9))
grid_plot_S$heights <- unit(rep(grid_row_height, 10), "cm")

# Save the grid as an image (e.g., PNG) with high resolution
ggsave("LME_grid_image_season.png", plot = grid_plot_S, width = 60, height = 100, units = "in", dpi = 300, limitsize = FALSE)


#Interaction Effects
#put in grid 
plot_list_IE <- list()
# Ensure Land_use and Season are factors
metadata$Land_use <- factor(metadata$Land_use, levels = c("Agri_ditch", "Mixed", "Forested"))
metadata$Season <- factor(metadata$Season, levels = c("Autumn", "Spring", "Summer"))
#all the pre-calculated pairwise pvalue from LME
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_AMIA_AMN_IE <- c("0.8130",
                                   "1",
                                   "1",
                                   "1",
                                   "1",
                                   "1",
                                   "0.0145",
                                   "1",
                                   "1")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_DOC_IE <- c("0.1501",
                              "0.0068",
                              "0.1071",
                              "0.0248",
                              "0.0093",
                              "0.0265",
                              "0.0065",
                              "0.0379",
                              "0.0424")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_DOC_MDL_IE <- c("0.0017",
                                  "0.0054",
                                  "0.0128",
                                  "1.8621e-06",
                                  "0.0003",
                                  "0.0003",
                                  "9.1988e-05",
                                  "0.0067",
                                  "0.0019")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_nitrate_IE <- c("0.0516",
                                  "1",
                                  "0.0291",
                                  "0.0503",
                                  "1",
                                  "0.1424",
                                  "0.3744",
                                  "1",
                                  "0.3215")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_TOC_IE <- c("0.1894",
                              "0.0070",
                              "0.1221",
                              "0.0268",
                              "0.0088",
                              "0.0272",
                              "0.0072",
                              "0.0378",
                              "0.0428")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_TOC_MDL_IE <- c("0.0017",
                                  "0.0054",
                                  "0.0128",
                                  "1.8621e-06",
                                  "0.0003",
                                  "0.0003",
                                  "9.1988e-05",
                                  "0.0067",
                                  "0.0019")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_TOTKN_IE <- c("0.0591",
                                "0.4546",
                                "0.3856",
                                "0.1662",
                                "1",
                                "0.3212",
                                "0.0011",
                                "1",
                                "0.3134")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_TOTKN_MDL_IE <- c("0.0408",
                                    "0.2084",
                                    "0.2140",
                                    "0.0152",
                                    "0.6410",
                                    "0.1440",
                                    "0.0014",
                                    "0.7480",
                                    "0.1035")
# Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_TOTPHO_IE <- c("0.5848",
                                 "1",
                                 "1",
                                 "1",
                                 "1",
                                 "1",
                                 "0.0084",
                                 "1")
 # Precalculated pairwise comparison p-values (replace with your actual p-values)
pairwise_p_values_TOTPHO_MDL_IE <- c("0.0476",
                                     "1",
                                     "0.3159",
                                     "0.0854",
                                     "1",
                                     "0.3199",
                                     "0.0021",
                                     "1",
                                     "0.1614")
# Create boxplot with Land_use and Season interaction

plot_list_IE[["AMIA_AMN_IE"]] <- ggplot(metadata, aes(x = interaction(Land_use, Season), y = AMIA_AMN, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$AMIA_AMN, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_AMIA_AMN_IE[1])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$AMIA_AMN, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_AMIA_AMN_IE[2])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$AMIA_AMN, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_AMIA_AMN_IE[3])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(metadata$AMIA_AMN, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_AMIA_AMN_IE[4])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(metadata$AMIA_AMN, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_AMIA_AMN_IE[5])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(metadata$AMIA_AMN, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_AMIA_AMN_IE[6])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(metadata$AMIA_AMN, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_AMIA_AMN_IE[7])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(metadata$AMIA_AMN, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_AMIA_AMN_IE[8])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(metadata$AMIA_AMN, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_AMIA_AMN_IE[9])),
              color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "AMIA_AMN") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#DOC
# Create boxplot with Land_use and Season interaction
plot_list_IE[["DOC_IE"]] <- ggplot(metadata, aes(x = interaction(Land_use, Season), y = DOC, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$DOC, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_DOC_IE[1])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$DOC, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_DOC_IE[2])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$DOC, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_DOC_IE[3])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(metadata$DOC, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_DOC_IE[4])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(metadata$DOC, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_DOC_IE[5])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(metadata$DOC, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_DOC_IE[6])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(metadata$DOC, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_DOC_IE[7])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(metadata$DOC, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_DOC_IE[8])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(metadata$DOC, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_DOC_IE[9])),
              color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "DOC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm"))+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#DOC_MDL

# Create boxplot with Land_use and Season interaction
plot_list_IE[["DOC_MDL_IE"]] <- ggplot(metadata, aes(x = interaction(Land_use, Season), y = DOC_MDL, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$DOC_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_DOC_MDL_IE[1])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$DOC_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_DOC_MDL_IE[2])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$DOC_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_DOC_MDL_IE[3])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(metadata$DOC_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_DOC_MDL_IE[4])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(metadata$DOC_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_DOC_MDL_IE[5])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(metadata$DOC_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_DOC_MDL_IE[6])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(metadata$DOC_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_DOC_MDL_IE[7])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(metadata$DOC_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_DOC_MDL_IE[8])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(metadata$DOC_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_DOC_MDL_IE[9])),
              color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "DOC_MDL") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#NITRATE
# Create boxplot with Land_use and Season interaction
plot_list_IE[["nitrate_IE"]] <- ggplot(metadata, aes(x = interaction(Land_use, Season), y = NITRATE, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$NITRATE, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_nitrate_IE[1])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$NITRATE, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_nitrate_IE[2])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$NITRATE, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_nitrate_IE[3])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(metadata$NITRATE, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_nitrate_IE[4])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(metadata$NITRATE, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_nitrate_IE[5])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(metadata$NITRATE, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_nitrate_IE[6])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(metadata$NITRATE, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_nitrate_IE[7])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(metadata$NITRATE, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_nitrate_IE[8])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(metadata$NITRATE, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_nitrate_IE[9])),
              color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "NITRATE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#TOC
# Create boxplot with Land_use and Season interaction
plot_list_IE[["TOC_IE"]] <- ggplot(metadata, aes(x = interaction(Land_use, Season), y = TOC, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$TOC, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOC_IE[1])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$TOC, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOC_IE[2])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$TOC, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOC_IE[3])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(metadata$TOC, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOC_IE[4])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(metadata$TOC, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOC_IE[5])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(metadata$TOC, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOC_IE[6])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(metadata$TOC, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOC_IE[7])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(metadata$TOC, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOC_IE[8])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(metadata$TOC, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOC_IE[9])),
              color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "TOC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#TOC_MDL

# Create boxplot with Land_use and Season interaction
plot_list_IE[["TOC_MDL_IE"]] <- ggplot(metadata, aes(x = interaction(Land_use, Season), y = TOC_MDL, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$TOC_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOC_MDL_IE[1])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$TOC_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOC_MDL_IE[2])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$TOC_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOC_MDL_IE[3])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(metadata$TOC_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOC_MDL_IE[4])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(metadata$TOC_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOC_MDL_IE[5])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(metadata$TOC_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOC_MDL_IE[6])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(metadata$TOC_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOC_MDL_IE[7])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(metadata$TOC_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOC_MDL_IE[8])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(metadata$TOC_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOC_MDL_IE[9])),
              color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "TOC_MDL") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#TOTKN

# Create boxplot with Land_use and Season interaction
plot_list_IE[["TOTKN_IE"]] <- ggplot(metadata, aes(x = interaction(Land_use, Season), y = TOTKN, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$TOTKN, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTKN_IE[1])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$TOTKN, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTKN_IE[2])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$TOTKN, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTKN_IE[3])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(metadata$TOTKN, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTKN_IE[4])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(metadata$TOTKN, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTKN_IE[5])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(metadata$TOTKN, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTKN_IE[6])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(metadata$TOTKN, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTKN_IE[7])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(metadata$TOTKN, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTKN_IE[8])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(metadata$TOTKN, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTKN_IE[9])),
              color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "TOTKN") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))



#TOTKN_MDL
# Create boxplot with Land_use and Season interaction
plot_list_IE[["TOTKN_MDL_IE"]] <- ggplot(metadata, aes(x = interaction(Land_use, Season), y = TOTKN_MDL, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$TOTKN_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTKN_MDL_IE[1])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$TOTKN_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTKN_MDL_IE[2])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$TOTKN_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTKN_MDL_IE[3])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(metadata$TOTKN_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTKN_MDL_IE[4])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(metadata$TOTKN_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTKN_MDL_IE[5])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(metadata$TOTKN_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTKN_MDL_IE[6])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(metadata$TOTKN_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTKN_MDL_IE[7])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(metadata$TOTKN_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTKN_MDL_IE[8])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(metadata$TOTKN_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTKN_MDL_IE[9])),
              color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "TOTKN_MDL") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))


#TOTPHO
# Create boxplot with Land_use and Season interaction
plot_list_IE[["TOTPHO_IE"]] <- ggplot(metadata, aes(x = interaction(Land_use, Season), y = TOTPHO, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$TOTPHO, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTPHO_IE[1])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$TOTPHO, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTPHO_IE[2])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$TOTPHO, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTPHO_IE[3])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(metadata$TOTPHO, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTPHO_IE[4])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(metadata$TOTPHO, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTPHO_IE[5])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(metadata$TOTPHO, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTPHO_IE[6])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(metadata$TOTPHO, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTPHO_IE[7])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(metadata$TOTPHO, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTPHO_IE[8])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(metadata$TOTPHO, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTPHO_IE[9])),
              color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "TOTPHO") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#TOTPHO_MDL
# Create boxplot with Land_use and Season interaction
plot_list_IE[["TOTPHO_MDL_IE"]] <- ggplot(metadata, aes(x = interaction(Land_use, Season), y = TOTPHO_MDL, fill = Land_use)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_text(aes(x = 1, y = max(metadata$TOTPHO_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTPHO_MDL_IE[1])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 2, y = max(metadata$TOTPHO_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTPHO_MDL_IE[2])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 3, y = max(metadata$TOTPHO_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTPHO_MDL_IE[3])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 4, y = max(metadata$TOTPHO_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTPHO_MDL_IE[4])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 5, y = max(metadata$TOTPHO_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTPHO_MDL_IE[5])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 6, y = max(metadata$TOTPHO_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTPHO_MDL_IE[6])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 7, y = max(metadata$TOTPHO_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTPHO_MDL_IE[7])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 8, y = max(metadata$TOTPHO_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTPHO_MDL_IE[8])),
              color = "black", size = 6, vjust = -0.5) +
  geom_text(aes(x = 9, y = max(metadata$TOTPHO_MDL, na.rm = TRUE) + 1, label = paste("p =", pairwise_p_values_TOTPHO_MDL_IE[9])),
              color = "black", size = 6, vjust = -0.5) +
  scale_fill_manual(values = land_use_colors) +
  labs(x = "Land_use x Season", y = "TOTPHO_MDL") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

# Arrange the figures in a 2 by 5 grid
grid_row_height <- 20  # Adjust as needed

# Arrange all figures in a 5 by 2 grid with specified height for each row
grid_plot_IE <- do.call(grid.arrange, c(grobs = plot_list_IE, ncol = 2, nrow = 5))
grid_plot_IE$heights <- unit(rep(grid_row_height, 10), "cm")

# Save the grid as an image (e.g., PNG) with high resolution
ggsave("LME_grid_image_IE.png", plot = grid_plot_IE, width = 30, height = 100, units = "in", dpi = 300, limitsize = FALSE)



