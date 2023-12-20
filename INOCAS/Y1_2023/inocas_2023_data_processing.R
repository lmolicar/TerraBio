#
# inocas_2023_data_processing
# ***************************
#
# Purpose:
#          This analysis is based on soil samples collected at INOCAS and proc-
#          essed by EcoMol.
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
source("../../../../../RCode/R_Scripts/triplet_fixer.R") # adjusted for local access

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
# KD 2023-09-19: This will always be true--EcoMol changes column names (schema)
# each time they send us data. lookupColnames should include EM colnames
# (changes), TB colnames (constant), and keep (Y/N)
lookupColnames <- read.csv("lookupColnames.csv")
lookupSitenames <- read.csv("lookupSitenames.csv")

# KD 2023-09-19: EcoMol changes these each time. Thus this needs to change each
# time which is why it is in the "Ingest" section... 
# LM todo: Change infoColnames below to match the INOCASSamples2023 file...

# columns differ from this version
infoColnames <- c("N", "storage", "volumeSampleID",
                  "replicate", "replicateVolume",
                  "EcoMolID", "concentrationDNA_nguL",
                  "purityDNA")

# 2023 data
inocas2023Raw <- read.csv("INOCAS_2023_bioinfo_results.csv",
                        stringsAsFactors = F,
                        sep=",",
                        col.names = lookupColnames$TB_ColName)

# LM 2023-09-25: These names apparently do not match the content of the file.
# The latter contains the following data:
# N - Ok                   
# Sample.ID.CIAT - No
# Site - No
# Amplification.Success - No
# Replicate - Ok            
# ID.EcoMol. - Ok         
# X.ng.μL. - Ok
# A260.280 - Ok

# KD 2023-09-19: Bad import for inocas2023Info because infoColnames above is
# incorrect. Remember, EcoMol changes colnames for the purity etc. information
# each time. 

inocas2023Info <- read.csv("INOCASSamples2023.csv",
                      stringsAsFactors = F, 
                      col.names = infoColnames,
                      header = T)


# KD 2023-09-19: Why do these functions have no arguments? Also do not recommend
# this approach to creating tables or doing it like this within a function.
# lulut should be created outside the function and passed as an argument.
# Otherwise you have to change the function for each different dataset...
# in additon, some tables are not passed to this function, so this function will not
# work. The function has no idea what inocas2023Info is.

qc_plot_purity <- function(){
  
  lulut <- data.frame(rbind(
    c("CF", "Counterfactual"),
    c("F",  "Floresta"),
    c("I",  "Intervention")
  ))
  
  colnames(lulut) <- c("code", "lutype")
  
  # Purity was imported as character -> convert to numeric
  inocas2023Info <- inocas2023Info%>% mutate(purityDNA = as.numeric(str_trim(purityDNA)))
  
  # KD 2023-09-19: This was imported as character due to white space in the input table. 
  
  
  inocas2023Info %>% select(EcoMolID, purityDNA) %>% 
    mutate(code = str_split(EcoMolID, "_", simplify = T)[,4]) %>%
    left_join(., y = lulut, by = "code") %>%
    ggplot(aes(x = lutype, y = purityDNA)) + geom_boxplot() + 
    geom_jitter() + 
    labs(x = "Land use type", y = "Purity (nm)")
}

# KD 2023-09-19: Same comments as above.
qc_plot_concentration <- function(){
  
  lulut <- data.frame(rbind(
    c("CF", "Counterfactual"),
    c("F",  "Floresta"),
    c("I",  "Intervention")
  ))
  
  colnames(lulut) <- c("code", "lutype")
  
  # Concentration was imported as character -> convert to numeric
  inocas2023Info <- inocas2023Info%>% mutate(concentrationDNA_nguL = as.numeric(str_trim(concentrationDNA_nguL)))
  
  inocas2023Info %>% select(EcoMolID, concentrationDNA_nguL) %>% 
    mutate(code = str_split(EcoMolID, "_", simplify = T)[,4]) %>%
    left_join(., y = lulut, by = "code") %>%
    ggplot(aes(x = lutype, y = concentrationDNA_nguL)) + geom_boxplot() + 
    labs(x = "Land use type", y = expression("Concentration" ~ ng~mu~L^{-1}))
  
}

# KD 2023-09-19: the summarize isn't necessary unless there is filtering already done. it should recreate the column "sampleTotalAbd"...

qc_plot_sampleTotalAbundance <- function(){
  
  dstdir <- "D:/OneDrive - CGIAR/Documents/CALPSE/Terrabio/ABF_DEALS/INOCAS/BIODIVERSITY_INDICATORS/ANALYSIS"
  dstfile <- "sampleTotalAbundance_landuse_unfiltered.png"
  dstfile <- paste(dstdir, dstfile, sep = "/")
  inocas2023Raw %>%
    group_by(sample) %>%
    summarise(abundance = sum(asvAbsoluteAbundance)) %>%
    mutate(code = str_split(sample, "-", simplify = T)[,3]) %>%
    left_join(., y = lulut, by = "code") %>%
    ggplot(aes(x = lutype, y = abundance)) + geom_point() +
    labs(x = "Land use type", y = "Sample total abundance",
         caption = "Unfiltered data")
  ggsave(dstfile, width = 6, height = 4, units = "in")
  
}



## ----- Data cleaning and setup | 2023 -------------------------------

# Remove columns we don't need for analysis.
head(inocas2023Raw)

# KD 2023-09-19: This should be saved as the 'filtered' table. See
# horta_2023_data_processing for example. Also want to strip the '>' from the
# ASV header for downstream analysis. 

inocas2023Raw <- inocas2023Raw[ , lookupColnames$Keep == "Y"]
head(inocas2023Raw)
# Initially 44948 items


# Remove sites not meeting minimum library size
removedSites <- unique(inocas2023Raw$sample[inocas2023Raw$sampleTotalAbd <= minlibrarySize])
# 2023-09-26 - Result: removed sites = 0; 
#              all reached minimum library size

inocas2023Raw <- inocas2023Raw[ !(inocas2023Raw$sample %in% removedSites) , ]

stopifnot(length(unique(inocas2023Raw$sample[inocas2023Raw$sampleTotalAbd <= minlibrarySize])) == 0)

print(removedSites)
remove(removedSites)


# Remove ASVs that don't meet criteria.

inocas2023Raw <- inocas2023Raw[ inocas2023Raw$primerExpectedLength == "in range", ]
# 2023-09-26 - Result: it dropped 6717 lines; now 32461
# KD 2023-09-19: I see this as 37148; the table() supports this as 7800 FALSE and 37148 TRUE.

# 2023-09-26 - Results: It dropped 29482; now 2979, but all the sites
#              reported arthropods (we have data for th 23 samples)
inocas2023Raw <- inocas2023Raw[ inocas2023Raw$phylumBLASTn %in% phylum, ]
# KD 2023-09-19: I see this as 3388; table(inocas2023Raw$phylumBLASTn) supports this. The presentation has the correct number.


# KD 2023-09-19: why is this a function? multiple arguments are not passed to
# this function, so this function will not work. The function has no idea what
# inocas2023Raw or lulut are, for example...

qc_plot_asvAbsAbundance <- function(){
  
  inocas2023Raw %>% select(metadata_3, asvAbsoluteAbundance) %>%
    left_join(., y = lulut, by = c("metadata_3" = "code")) %>%
    ggplot(aes(x = lutype, y = asvAbsoluteAbundance)) +
    geom_boxplot() + labs(x = "Land use type",
                          y = "ASV Absolute Abundance",
                          caption = "Filtered data")
  
}

# KD 2023-09-19: Same comment as above.

qc_plot_sampleTotalAbundance <- function(){
  
  inocas2023Raw %>%
    group_by(sample) %>%
    summarise(abundance = sum(asvAbsoluteAbundance)) %>%
    mutate(code = str_split(sample, "-", simplify = T)[,3]) %>%
    left_join(., y = lulut, by = "code") %>%
    ggplot(aes(x = lutype, y = abundance)) + geom_point() +
    labs(x = "Land use type", y = "Sample total abundance",
         caption = "Filtered data")
  
}

# KD 2023-09-19: Same comment as above.

qc_plot_counts_landuse <- function(){
  
  inocas2023Raw %>% select(sample, classBLASTn) %>%
    mutate(code = str_split(sample, "-", simplify = TRUE)[,3]) %>%
    left_join(., y = lulut, by = "code") %>%
    filter(lutype == "Counterfactual") %>%
    ggplot(aes(x = classBLASTn)) + geom_bar() + coord_flip() + 
    labs(x = "BLASTn Class", title = "Counterfactual")
  
}

# KD 2023-09-19: Same comment as above, but with the addition of "curlanduse"...

qc_plot_asvcounts_class_landuse <- function(){
  
  curlanduse <- "Intervention"
  dstdir <- "D:/OneDrive - CGIAR/Documents/CALPSE/Terrabio/ABF_DEALS/INOCAS/BIODIVERSITY_INDICATORS/ANALYSIS"
  dstfile <- paste0("blastn-class_asvabundance_",
                    curlanduse, ".png")
  dstfile <- paste(dstdir, dstfile, sep = "/")
  inocas2023Raw %>% select(sample, classBLASTn, asvAbsoluteAbundance) %>%
    mutate(code = str_split(sample, "-", simplify = TRUE)[,3]) %>%
    left_join(., y = lulut, by = "code") %>%
    filter(lutype == curlanduse) %>%
    group_by(classBLASTn) %>%
    summarize(asvCount = sum(asvAbsoluteAbundance)) %>%
    ggplot(aes(x = classBLASTn, y = asvCount)) +
    geom_bar(stat = "identity") +
    coord_flip() + 
    labs(x = "BLASTn Class",
         y = "ASV Absolute Abundance",
         title = curlanduse)
  
  ggsave(dstfile, units = "in", width = 6, height = 3)
  
}


# KD 2023-09-19: Same comment as above, but with the addition of "curlanduse"...


qc_plot_asvcounts_order_landuse <- function(){
  
  curlanduse <- "Floresta"
  dstdir <- "D:/OneDrive - CGIAR/Documents/CALPSE/Terrabio/ABF_DEALS/INOCAS/BIODIVERSITY_INDICATORS/ANALYSIS"
  dstfile <- paste0("blastn-order_asvabundance_",
                    curlanduse, ".png")
  dstfile <- paste(dstdir, dstfile, sep = "/")
  inocas2023Raw %>% select(sample, order_BLASTn, asvAbsoluteAbundance) %>%
    mutate(code = str_split(sample, "-", simplify = TRUE)[,3]) %>%
    left_join(., y = lulut, by = "code") %>%
    filter(lutype == curlanduse) %>%
    group_by(order_BLASTn) %>%
    summarize(asvCount = sum(asvAbsoluteAbundance)) %>%
    ggplot(aes(x = order_BLASTn, y = asvCount)) +
    geom_bar(stat = "identity") +
    coord_flip() + 
    labs(x = "BLASTn Order",
         y = "ASV Absolute Abundance",
         title = curlanduse)
  
  ggsave(dstfile, units = "in", width = 5, height = 6)
  
}


# Luis: Checks in Apui 2022 case - Should apply in INOCAS? Always?

# KD 2023-09-19: Yes, these checks should always be completed. This is an
# important part of data filtering, with asvAbsoluteAbundance being the prefered
# metric based on what we discussed in the previous learning session. See also
# horta_2023_data_processing. Also suggest remove the "if (TRUE)".

if (TRUE){
  
  plot(
    inocas2023Raw$asvAbsoluteAbundance,
    inocas2023Raw$relativeAbundanceOnSample,
    ylim = c(0, 0.1),
    xlim = c(0, 5000),
    xlab = "ASV absolute abundance",
    ylab = "Relative abundance on sample"
  )
  abline(h = minRelativeAbund, col = "red")
  abline(v = minAbsoluteAbund, col = "red")
  
  print("Number of ASV under min relative abundance:")
  print(sum(inocas2023Raw$relativeAbundanceOnSample < minRelativeAbund))
  print("Number of ASV under min absolute abundance:")
  print(sum(inocas2023Raw$asvAbsoluteAbundance < minAbsoluteAbund))
  print(sum(inocas2023Raw$asvAbsoluteAbundance > minAbsoluteAbund))
  
  inocas2023Raw <- inocas2023Raw[ inocas2023Raw$asvAbsoluteAbundance > minAbsoluteAbund, ]
  
  # KD 2023-09-19: 2702 rows following filtering
  
  # KD 2023-09-19: This change is required b/c this code doesn't use inocas2023Filt. 
  
  inocasASV <- inocas2023Raw
  
  # What is this for? Where is it used?
  # Species abundance (absolute) per sample
  
  # KD 2023-09-19: This is one of the two main site-species matrices used for
  # the indicator calculation... should be calculated after all filtering steps
  # are complete.
  
  inocasMatrix <- ez.matrify(inocasASV,
                             species.name = "ASVHeader",
                             site.name = "sample",
                             abundance = "asvAbsoluteAbundance")
  
}


# Create a dataset where the land use type are combined by replicate
# In Horta 2022 data:
# metadata_5 = preservation {"b", "s"}
# primerName = primer {"R1", "R2"}
# metadata_4 = treatment {"CF", "Co", "Re", "Sy"}
# Thus, sample letter was something like: b-R1-CF
# hortaASV$sampleLetter <- paste0(hortaASV$metadata_5,
#                                 "-",
#                                 hortaASV$primerName,
#                                 "-",
#                                 str_sub(hortaASV$metadata_4, 1, 2))
#
# For INOCAS it'simpler: we have only one preservation method, and one
# primer. We are computing the absolute abundance per land use, i.e.,
# per treatment, which is stored in "metadata_3"
#
inocasASV$sampleLetter <- inocasASV$metadata_3
inocasLetter <- inocasASV %>%
  dplyr::select(sampleLetter, ASVHeader, asvAbsoluteAbundance) %>%
  group_by(sampleLetter, ASVHeader) %>%
  summarise(abundance = sum(asvAbsoluteAbundance))

inocasMatrixLetter <- ez.matrify(inocasLetter, species.name = "ASVHeader",
                                site.name = "sampleLetter", abundance = "abundance")

#test to make sure everything got in
any((colSums(inocasMatrixLetter)-colSums(inocasMatrix)) > 0 )

# KD 2023-09-19: Recommend removing data tables note needed for the indicator calculations because it reduces memory load etc.

#remove(inocas2023Raw)


# KD 2023-09-19: These were to-dos for previous code, recommend removing.

# to do list:
#   + add column names for 2023
#   + import 2023 data (even if wrong) to build out
#   + need to look at which ASVs should be combined or not...






