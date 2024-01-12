#
# inocas_2023_data_processing
# ***************************
#
# Purpose:
#          This analysis is based on soil samples collected at INOCAS and proc-
#          essed by EcoMol in 2023.
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

# LM 2023-09-25: Check column names and those in the lookup table (LUT);
#                there were more columns in the LUT than in the data.
#                I adjusted the LUT to reflect the columns in data.
#                This may have to be adjusted depending on whether the
#                new version of the data includes columns that are
#                required.
#                Original version created by Karen was stored as
#                lookupColnames_v1.csv.bak
#
# KD 2023-09-19: This will always be true--EcoMol changes column names (schema)
# each time they send us data. lookupColnames should include EM colnames
# (changes), TB colnames (constant), and keep (Y/N)

lookupColnames <- read.csv("lookupColnames.csv")
lookupSitenames <- read.csv("lookupSitenames.csv")

# KD 2023-09-19: EcoMol changes these each time. Thus this needs to change each
# time which is why it is in the "Ingest" section... 

# columns differ from this version
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


##---- Quality control checks

# KD 2023-09-19: Why do these functions have no arguments? Also do not recommend
# this approach to creating tables or doing it like this within a function.
# lulut should be created outside the function and passed as an argument.
# Otherwise you have to change the function for each different dataset...
# in additon, some tables are not passed to this function, so this function will not
# work. The function has no idea what inocas2023Info is.

#-- EcoMol ID contains a code for the treatment. We will associate that
# code with a nice name for the treatment.

lulut <- data.frame(code = c("CF", "F", "I"),
                    lutype = c("Counterfactual", "Reference", "Intervention"))

#-- Sample purity per treatment (land use)
  
# Purity was imported as character -> convert to numeric due to white space
# in the input table (detected by KD 2023-12-19)
inocas2023Info <- inocas2023Info %>%
  mutate(purityDNA = as.numeric(str_trim(purityDNA)))


inocas2023Info %>% select(EcoMolID, purityDNA) %>% 
  mutate(code = str_split(EcoMolID, "_", simplify = T)[,4]) %>%
  left_join(., y = lulut, by = "code") %>%
  ggplot(aes(x = lutype, y = purityDNA)) + geom_boxplot() + 
  geom_jitter() + 
  labs(x = "Land use type", y = "Purity (nm)")

#-- Sample concentration per land use
  
# Concentration was imported as character -> convert to numeric
inocas2023Info <- inocas2023Info%>% mutate(concentrationDNA_nguL = as.numeric(str_trim(concentrationDNA_nguL)))

inocas2023Info %>% select(EcoMolID, concentrationDNA_nguL) %>% 
  mutate(code = str_split(EcoMolID, "_", simplify = T)[,4]) %>%
  left_join(., y = lulut, by = "code") %>%
  ggplot(aes(x = lutype, y = concentrationDNA_nguL)) + geom_boxplot() + 
  labs(x = "Land use type", y = expression("Concentration" ~ ng~mu~L^{-1}))

#-- Sample total abundance per land use
  
dstdir <- "D:/OneDrive - CGIAR/Documents/CALPSE/Terrabio/ABF_DEALS/INOCAS/BIODIVERSITY_INDICATORS/ANALYSIS"
dstfile <- "sampleTotalAbundance_landuse_unfiltered.png"
dstfile <- paste(dstdir, dstfile, sep = "/")
inocas2023Raw %>%
  select(sample, sampleTotalAbd) %>% distinct() %>%
  mutate(code = str_split(sample, "-", simplify = T)[,3]) %>%
  left_join(., y = lulut, by = "code") %>%
  ggplot(aes(x = lutype, y = sampleTotalAbd)) + geom_point() +
  labs(x = "Land use type", y = "Sample total abundance",
       caption = "Unfiltered data")
ggsave(dstfile, width = 6, height = 4, units = "in")



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

#-- ASV absolute abundance per land use
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

#-- Sample total abundance following filtering
inocas2023Filtered %>%
  group_by(sample) %>%
  summarise(sampleTotalAbd = sum(asvAbsoluteAbundance)) %>%
  mutate(code = str_split(sample, "-", simplify = T)[,3]) %>%
  left_join(., y = lulut, by = "code") %>%
  ggplot(aes(x = lutype, y = sampleTotalAbd)) + geom_point() +
  labs(x = "Land use type", y = "Sample total abundance",
       caption = "Filtered data")


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
# the indicator calculation should be calculated after all filtering steps are complete (KD 2023-12-19).

# LM 2024-01-12: I realized that here this is actually for replicates. I recomputed this in inocas_2023_indicators, where I use an ID for each sample (site)
# inocasMatrix <- ez.matrify(inocas2023Filtered,
#                            species.name = "ASVHeader",
#                            site.name = "sample",
#                            abundance = "asvAbsoluteAbundance")


# For INOCAS it's simpler--we have only one preservation method, and one
# primer. We are computing a treatment-species matrix, i.e., the absolute
# abundance per land use (or treatment), which is stored in "metadata_3"
#

# LM 2024-01-11: I think the previous code (below) could be simplified. The function ez.matrify summarizes the data, right?

inocas2023Filtered$sampleLetter <- inocas2023Filtered$metadata_3
# inocasLetter <- inocas2023Filtered %>%
#   dplyr::select(sampleLetter, ASVHeader, asvAbsoluteAbundance) %>%
#   group_by(sampleLetter, ASVHeader) %>%
#   summarise(abundance = sum(asvAbsoluteAbundance))

inocasMatrixLetter <- ez.matrify(inocas2023Filtered,
                                 species.name = "ASVHeader",
                                 site.name = "sampleLetter",
                                 abundance = "asvAbsoluteAbundance")

#test to make sure everything got in
any((colSums(inocasMatrixLetter)-colSums(inocasMatrix)) > 0 )

# Remove raw data to free up resources
remove(inocas2023Raw)
