#
# inocas_2023_data_processing
# ***************************
#
# Purpose:
#          analysis of soil samples collected at INOCAS and processed by EcoMol
#          in 2023.
#          Takes the eDNA species table provided by EcoMol and prepares them for
#          analysis.
#
#
# Date             Written by                 Description of changes
# ***********  *****************  **********************************************
# Sep 2023     Karen Dyson        Original code.
#
# Dec 14 2023  Luis Molina        Adapted to read INOCAS data 2023.
#
# Dec 19 2023  Karen Dyson        Comments added/code review
#
# 


## ----- Data ingestion ---------------------------------------

# libraries
library(ggplot2)
library(tidyr)
library(vegan)
library(stringr)
library(compositions)
library(zCompositions)
library(dplyr)

# Ingest codes
source("../../allianceBranding.R")
source("../../functions.R")
source("../../multiyear_functions.R")
source("../../../RCode/KDyson_R_Scripts/triplet_fixer.R") # adjusted for local access

# Special flags
saveplots <- FALSE

# script variable definitions
minlibrarySize = 5000
minRelativeAbund = 0.05
minAbsoluteAbund = 5
rare = 50

#Reitmeier, S., Hitch, T.C., Treichel, N., Fikas, N., Hausmann, B., Ramer-Tait,
#A.E., Neuhaus, K., Berry, D., Haller, D., Lagkouvardos, I. and Clavel, T.,
#2021. Handling of spurious sequences affects the outcome of high-throughput 16S
#rRNA gene amplicon profiling. ISME Communications, 1(1), pp.1-12.


phylum = c("Arthropoda")
#phylum = c("Annelida", "Nematoda", "Platyhelminthes", "Arthropoda", "Mollusca") #worms and insects



## First ingest 2023 data.

## common data

# Check column names and those in the lookup table (LUT); this will always be
# true--EcoMol changes column names (schema) each time they send data. The table
# lookupColnames should include:
# - EM colnames (it changes)
# - TB colnames (constant) and
# - keep (Y/N)

lookupColnames <- read.csv("lookupColnames.csv")
lookupSitenames <- read.csv("lookupSitenames.csv")

# KD 2023-09-19: EcoMol changes these each time. Thus this needs to change each
# time which is why it is in the "Ingest" section... 
infoColnames <- c("N", "ciatID", "treatment",
                  "ampliSuccess", "replNo",
                  "EcoMolID", "concentrationDNA_nguL",
                  "purityDNA")

# 2023 data
inocas2023Raw <- read.csv("INOCAS_2023_bioinfo_results.csv",
                        stringsAsFactors = F,
                        sep=",",
                        col.names = lookupColnames$TB_ColName)

# Always check column names (KD 2023-09-19): EcoMol changes colnames for the purity etc. information each time.
inocas2023Info <- read.csv("INOCASSamples2023.csv",
                      stringsAsFactors = F, 
                      col.names = infoColnames,
                      header = T)


#----  Quality control checks ----

## EcoMol ID contains a code for the treatment. We will associate that
# code with a nice name for the treatment.

lulut <- data.frame(code = c("CF", "F", "I"),
                    lutype = c("Counterfactual", "Reference", "Intervention"))

#-- Sample purity per treatment (land use)
  
# Purity was imported as character -> convert to numeric due to white space
# in the input table (detected by KD 2023-12-19)
inocas2023Info <- inocas2023Info %>%
  mutate(purityDNA = as.numeric(str_trim(purityDNA)))

# Checking samples that arrived in the lab without preservation method.
#Samples identified manually. In the future we can include a column whether sam-
# ples had a preservation method.
samplesIDNoPreserve <- c("EM135c2_Tbio3_3_I_R1",
                         "EM135c2_Tbio3_3_F_R1",
                         "EM135c2_Tbio3_3_CF_R2",
                         "EM135c2_Tbio3_3_CF_R3")


inocas2023Info %>% select(EcoMolID, purityDNA) %>% 
  mutate(code = str_split(EcoMolID, "_", simplify = T)[,4]) %>%
  left_join(., y = lulut, by = "code") %>%
  mutate(preserve = if_else(EcoMolID %in% samplesIDNoPreserve, "No preservation", "Silica")) %>% filter(preserve == "") %>%
  ggplot(aes(x = lutype, y = purityDNA)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(col=lutype),width = 0.1) + 
  geom_text(aes(label = preserve), nudge_x = 0.2, size = 5) +
  labs(color = "Site Type", y = "Purity (nm)",
       x = "Site Type")+
  scale_color_manual(values = supportingColorPalette[c(1,3,4)])+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14))

inocas2023Info %>% select(EcoMolID, purityDNA) %>% 
  mutate(code = str_split(EcoMolID, "_", simplify = T)[,4]) %>%
  left_join(., y = lulut, by = "code") %>%
  mutate(preserve = if_else(EcoMolID %in% samplesIDNoPreserve, "No preservation", "Silica")) %>%
  ggplot(aes(x = lutype, y = purityDNA)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(col=lutype, shape = preserve), width = 0.3, size = 4) +
  labs(color = "Site Type", y = "Purity (nm)",
       x = "Site Type")+
  scale_color_manual(values = supportingColorPalette[c(1,3,4)])+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14))


if (saveplots)
  ggsave("purity_per_landuse.png", width = 8, height = 5, units = "in", dpi = 300)


#-- Sample concentration per land use
  
# Concentration was imported as character -> convert to numeric
inocas2023Info <- inocas2023Info%>% mutate(concentrationDNA_nguL = as.numeric(str_trim(concentrationDNA_nguL)))

inocas2023Info %>% select(EcoMolID, concentrationDNA_nguL) %>% 
  mutate(code = str_split(EcoMolID, "_", simplify = T)[,4]) %>%
  left_join(., y = lulut, by = "code") %>%
  mutate(preserve = if_else(EcoMolID %in% samplesIDNoPreserve, "No preservation", "Silica")) %>%
  ggplot(aes(x = lutype, y = concentrationDNA_nguL)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(col=lutype, shape = preserve), width = 0.3, size =4) + 
  labs(color = "Land use type", x = "Land use type", y = expression("Concentration" ~ ng~mu~L^{-1})) +scale_color_manual(values = supportingColorPalette[c(1,3,4)])+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14))

if (saveplots)
  ggsave("concentration_per_landuse.png", width = 8, height = 5, units = "in", dpi = 300)

#-- Sample total abundance per land use
  
dstdir <- "D:/OneDrive - CGIAR/Documents/CALPSE/Terrabio/ABF_DEALS/INOCAS/BIODIVERSITY_INDICATORS/ANALYSIS"
dstfile <- "sampleTotalAbundance_landuse_unfiltered.png"
#dstfile <- paste(dstdir, dstfile, sep = "/")
inocas2023Raw %>%
  select(sample, sampleTotalAbd) %>% distinct() %>%
  mutate(code = str_split(sample, "-", simplify = T)[,3]) %>%
  left_join(., y = lulut, by = "code") %>%
  ggplot(aes(x = lutype, y = sampleTotalAbd)) + geom_point() +
  labs(x = "Land use type", y = "Sample total abundance",
       caption = "Unfiltered data")

if (saveplots)
  ggsave(dstfile, width = 6, height = 4, units = "in")

#--- sample abundance per land use (differences in preservation method)

inocas2023Raw %>%
  select(sample, sampleTotalAbd) %>% 
  mutate(preserve = if_else(sample %in% str_replace(str_replace_all(samplesIDNoPreserve, "_", "-"), "-Tbio3-", "-"), "No preservation", "Silica")) %>%
  distinct() %>%
  mutate(code = str_split(sample, "-", simplify = T)[,3]) %>%
  left_join(., y = lulut, by = "code") %>%
  ggplot(aes(x = lutype, y = sampleTotalAbd)) + geom_point(aes(col = preserve), size = 3) +
  labs(x = "Land use type", y = "Number of reads",
       caption = "Unfiltered data", col = "Preservation")+
  ggtitle("Sample absolute abundance")+
  theme(plot.title = element_text(hjust = 0.5))

if (saveplots)
  ggsave(dstfile, width = 6, height = 4, units = "in")


#-- ASV absolute abundance per land use (unfiltered data)
inocas2023Raw %>% select(metadata_3, asvAbsoluteAbundance) %>%
  left_join(., y = lulut, by = c("metadata_3" = "code")) %>%
  ggplot(aes(x = lutype, y = asvAbsoluteAbundance)) +
  geom_boxplot() + labs(x = "Land use type",
                        y = "Abundance",
                        caption = "Unfiltered data")+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14)) + 
  ggtitle("ASV absolute abundance") + 
  theme(plot.title = element_text(hjust = 0.5))

if (saveplots)
  ggsave("asvAbsoluteAbundance_landuse_unfiltered.png", width = 6,
       height = 4, units = "in")

# What's the composition of raw sample data (i.e., unfiltered)
# We know we have man unidentified ASVs -> remove them first.

#***** Order level
asvCountUnfiltered <- inocas2023Raw %>% select(phylumBLASTn, order_BLASTn) %>%
  #mutate(code = str_split(sample, "-", simplify = TRUE)[,3]) %>%
  filter(phylumBLASTn == "Arthropoda") %>% select(order_BLASTn) %>% table() %>% as.data.frame() %>% arrange(desc(Freq))

asvCountUnfiltered %>% ggplot(aes(reorder(order_BLASTn, -Freq), Freq)) + geom_col() + coord_flip() + ggtitle("All samples") +
  labs(x = "Order",
       y = "ASV Count")+
  theme(plot.title = element_text(hjust = 0.5))

asvCountUnfilteredNoPreserve <- inocas2023Raw %>% select(phylumBLASTn, order_BLASTn, metadata_3, sample) %>%
  left_join(., y = lulut, by = c("metadata_3" = "code")) %>%
  mutate(preserve = if_else(sample %in% str_replace(str_replace_all(samplesIDNoPreserve, "_", "-"), "-Tbio3-", "-"), "No preservation", "Silica")) %>% 
  filter(preserve == "No preservation") %>%
  filter(phylumBLASTn == "Arthropoda") %>% select(order_BLASTn) %>% table() %>% as.data.frame() %>% arrange(desc(Freq))

asvCountUnfilteredNoPreserve <- asvCountUnfilteredNoPreserve %>% rename(flag = Freq) %>% mutate(flag = "No silica") %>% right_join(., asvCountUnfiltered) %>% mutate(flag = if_else(is.na(flag), "Missing", flag))

flag_sorted <- asvCountUnfilteredNoPreserve %>% select(Freq, flag) %>% arrange(-Freq)

#>>> trying to apply colors to individual ticks in y axis
#>>> This link was useful https://community.rstudio.com/t/axis-labels-with-individual-colors/77848
ticks_colors <- if_else(flag_sorted$flag == "Missing", "coral1", "#505050")

asvCountUnfilteredNoPreserve %>% ggplot(aes(reorder(order_BLASTn, -Freq), Freq, fill = flag)) + geom_col() + coord_flip() + ggtitle("All samples") +
  labs(x = "Order",
       y = "ASV Count")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(colour = ticks_colors)) +
  scale_fill_manual(values=c("coral1", "#505050"))


if (saveplots)
  ggsave("asv_count_order_unfiltered.png", width = 6, height = 6, units = "in", dpi = 300)


#-- Let's check ASV abundance of the samples that did not have silica
inocas2023Raw %>% select(metadata_3, asvAbsoluteAbundance, sample) %>%
  left_join(., y = lulut, by = c("metadata_3" = "code")) %>%
  mutate(preserve = if_else(sample %in% str_replace(str_replace_all(samplesIDNoPreserve, "_", "-"), "-Tbio3-", "-"), "No preservation", "Silica")) %>% 
  filter(preserve == "No preservation") %>%
  ggplot(aes(x = lutype, y = asvAbsoluteAbundance)) +
  geom_point() + labs(x = "Land use type",
                        y = "Abundance",
                        caption = "Unfiltered data")+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14)) + 
  ggtitle("ASV absolute abundance") + 
  theme(plot.title = element_text(hjust = 0.5))

if (saveplots)
  ggsave("asvAbsoluteAbundance_landuse_unfiltered_NoSilica.png", width = 6,
         height = 4, units = "in")

#***** And what's its composition?
# asvCountUnfilteredNoPreserve <- inocas2023Raw %>% select(phylumBLASTn, order_BLASTn, metadata_3, sample) %>%
#   left_join(., y = lulut, by = c("metadata_3" = "code")) %>%
#   mutate(preserve = if_else(sample %in% str_replace(str_replace_all(samplesIDNoPreserve, "_", "-"), "-Tbio3-", "-"), "No preservation", "Silica")) %>% 
#   filter(preserve == "No preservation") %>%
#   filter(phylumBLASTn == "Arthropoda") %>% select(order_BLASTn) %>% table() %>% as.data.frame() %>% arrange(desc(Freq))

asvCountUnfilteredNoPreserve %>% ggplot(aes(reorder(order_BLASTn, -Freq), Freq)) + geom_col(fill = "#505050") + coord_flip() + ggtitle("Samples without silica") +
  labs(x = "Order",
       y = "ASV Count")+
  theme(plot.title = element_text(hjust = 0.5))


if (saveplots)
  ggsave("asv_count_order_unfiltered_NoSilica.png", width = 6, height = 6, units = "in", dpi = 300)


## ----- Data cleaning and setup | 2023 -------------------------------

# Remove columns we don't need for analysis.
head(inocas2023Raw)

# Save a filtered version of the input table
inocas2023Filtered <- inocas2023Raw[ , lookupColnames$Keep == "Y"]

# Remove the '>' from the ASV header for donwstream analysis
inocas2023Filtered$ASVHeader <- str_sub(inocas2023Filtered$ASVHeader, 2)
head(inocas2023Filtered)
# Initially 44948 items


# Remove sites not meeting minimum library size
removedSites <- unique(inocas2023Filtered$sample[inocas2023Filtered$sampleTotalAbd <= minlibrarySize])
# 2023-09-26 - Result: removed sites = 0; 
#              all reached minimum library size

inocas2023Filtered <- inocas2023Filtered[ !(inocas2023Filtered$sample %in% removedSites) , ]

stopifnot(length(unique(inocas2023Filtered$sample[inocas2023Filtered$sampleTotalAbd <= minlibrarySize])) == 0)

print(removedSites)
remove(removedSites)

# Remove ASVs that don't meet criteria.

#-- Filter for primer length
inocas2023Filtered <- inocas2023Filtered[ inocas2023Filtered$primerExpectedLength == "in range", ]
# table() reports 7800 lines dropped; 37148 lines remaining

#-- Filter for target phylum(s)
inocas2023Filtered <- inocas2023Filtered[ inocas2023Filtered$phylumBLASTn %in% phylum, ]
# Result: 3388 lines remain after filtering


# Exploratory data analysis ----

#-- ASV absolute abundance per land use (filtered data)
inocas2023Filtered %>% select(metadata_3, asvAbsoluteAbundance) %>%
  left_join(., y = lulut, by = c("metadata_3" = "code")) %>%
  ggplot(aes(x = lutype, y = asvAbsoluteAbundance)) +
  geom_boxplot() + labs(x = "Land use type",
                        y = "Abundance",
                        caption = "Filtered data")+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14)) + 
  ggtitle("ASV absolute abundance") + 
  theme(plot.title = element_text(hjust = 0.5))

if (saveplots)
  ggsave("asvAbsoluteAbundance_landuse_filtered.png", width = 6,
       height = 4, units = "in")


# Some high values in intervention and reference sites - What are they?
# Values in question are those > 5000 in intervention and > 1000 in reference
t(inocas2023Raw %>% select(metadata_3,
                           asvAbsoluteAbundance,
                           colnames(.)[grepl("BLASTn", colnames(.))]) %>%
    left_join(., y = lulut, by = c("metadata_3" = "code")) %>%
    filter(lutype == "Intervention") %>% filter(asvAbsoluteAbundance > 5000))


#-- Sample total abundance following filtering
inocas2023Filtered %>%
  group_by(sample) %>%
  summarise(sampleTotalAbd = sum(asvAbsoluteAbundance)) %>%
  mutate(code = str_split(sample, "-", simplify = T)[,3]) %>%
  left_join(., y = lulut, by = "code") %>%
  ggplot(aes(x = lutype, y = sampleTotalAbd)) + geom_point() +
  labs(x = "Land use type", y = "Abundance",
       caption = "Filtered data")+
  ggtitle("Sample total abundance")+
  theme(plot.title = element_text(hjust = 0.5))


inocas2023Filtered %>% 
  mutate(preserve = if_else(sample %in% str_replace(str_replace_all(samplesIDNoPreserve, "_", "-"), "-Tbio3-", "-"), "No preservation", "Silica")) %>%
  group_by(sample, preserve) %>%
  summarise(sampleTotalAbd = sum(asvAbsoluteAbundance)) %>%
  mutate(code = str_split(sample, "-", simplify = T)[,3]) %>%
  left_join(., y = lulut, by = "code") %>%
  ggplot(aes(x = lutype, y = sampleTotalAbd, col = preserve)) + 
  geom_point(size = 3) +
  labs(x = "Land use type", y = "Abundance",
       caption = "Filtered data", col = "Preservation") +
  ggtitle("Sample total abundance")+
  theme(plot.title = element_text(hjust = 0.5))

if (saveplots)
  ggsave("sampleTotalAbundance_landuse_filtered.png", width = 6, height = 4, units = "in")

# Checks below should always be completed (Karen Dyson, 2023-12-19).
# Note that asvAbsoluteAbundance is the preferred metric (check learning
# session with Karen). Code following horta_2023_data_processing.

plot(
  inocas2023Filtered$asvAbsoluteAbundance,
  inocas2023Filtered$relativeAbundanceOnSample,
  #ylim = c(0, 0.05),
  #xlim = c(0, 50),
  xlab = "ASV absolute abundance",
  ylab = "Relative abundance on sample"
)

abline(v = minAbsoluteAbund, col = "red")

print("Number of filtered ASVs under minimum absolute abundance:")
print(sum(inocas2023Filtered$asvAbsoluteAbundance < minAbsoluteAbund))

inocas2023Filtered <- inocas2023Filtered[ inocas2023Filtered$asvAbsoluteAbundance > minAbsoluteAbund, ]

# Result: 2702 rows following filtering (KD 2023-09-19)

# ---- Site-species matrices ----

# This is one of the two main site-species matrices used for
# the indicator calculation should be calculated after all filtering steps are
# complete (KD 2023-12-19).

# LM 2024-01-12: I realized that here this is actually for replicates (samples).
# I recomputed this in inocas_2023_indicators, where I use an ID for each sample
# (site)
inocasMatrix <- ez.matrify(inocas2023Filtered,
                           species.name = "ASVHeader",
                           site.name = "sample",
                           abundance = "asvAbsoluteAbundance")


# For INOCAS it's simpler--we have only one preservation method, and one
# primer. We are computing a treatment-species matrix, i.e., the absolute
# abundance per land use (or treatment), which is stored in "metadata_3"
#

inocas2023Filtered$sampleLetter <- inocas2023Filtered$metadata_3

# LM 2024-01-11: The function ez.matrify summarizes the data. No need to
# summarize.
inocasMatrixLetter <- ez.matrify(inocas2023Filtered,
                                 species.name = "ASVHeader",
                                 site.name = "sampleLetter",
                                 abundance = "asvAbsoluteAbundance")

# Test to make sure everything got in
any((colSums(inocasMatrixLetter)-colSums(inocasMatrix)) > 0 )

# Remove raw data to free up resources
remove(inocas2023Raw)
