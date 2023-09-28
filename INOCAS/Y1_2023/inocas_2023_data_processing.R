# This analysis is based on soil samples collected at Horta and processed
# by EcoMol.

# This file take the eDNA species table provided by EcoMol and prepares them for
# analysis.

# Code written Sept. 2023 by Karen Dyson


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
source("../allianceBranding.R")
source("../functions.R")
source("../multiyear_functions.R")
source("../../R_Scripts/triplet_fixer.R") # adjusted for local access

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



## First ingest 2022 and 2023 data.

## common data

# LM 2023-09-25: Check column names and those in the lookup table (LUT);
#                there were more columns in the LUT than in the data.
#                I adjusted the LUT to reflect the columns in data.
#                This may have to be adjusted depending on whether the
#                new version of the data includes columns that are
#                required.
#                Original version created by Karen was stored as
#                lookupColnames_v1.csv.bak
lookupColnames <- read.csv("lookupColnames.csv")
lookupSitenames <- read.csv("lookupSitenames.csv")


# columns differ from this version
# infoColnames <- c("N", "storage", "volumeSampleID",
#                   "replicate", "replicateVolume",
#                   "EcoMolID", "concentrationDNA_nguL",
#                   "purityDNA")
infoColnames <- c(
  "N",
  "ciatID",
  "site",
  "ampSuccess",
  "preservation",
  "replicate",
  "ecomolID",
  "concentrationDNA_nguL",
  "purityDNA")

# 2023 data
inocas2023Raw <- read.csv("Y1_2023/INOCAS_2023_bioinfo_results.csv",
                        stringsAsFactors = F,sep=",",
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
inocas2023Info <- read.csv("Y1_2023/INOCASSamples2023.csv",
                      stringsAsFactors = F,
                      col.names = infoColnames,
                      header = T)




## ----- Data cleaning and setup | 2023 -------------------------------

# Remove columns we don't need for analysis.
head(inocas2023Raw)
inocas2023Raw <- inocas2023Raw[ , lookupColnames$Keep == "Y"]
head(inocas2023Raw)
# Initially 44948 items

# Check for preservation method - Remove samples that arrived
# at the lab without preservation
inocas2023Info$sample <- gsub("Tbio3", "EM135c2", x = gsub("_", "-", x = inocas2023Info$ciatID))
sitesNoPreservation <- inocas2023Info$sample[inocas2023Info$preservation == "none"]
inocas2023Raw <- inocas2023Raw[ ! inocas2023Raw$sample %in% sitesNoPreservation, ]
# it removed 5770 items (39170)

# Remove sites not meeting minimum library size
removedSites <- unique(inocas2023Raw$sample[inocas2023Raw$sampleTotalAbd <= minlibrarySize])
# 2023-09-26 - Result: removed sites = 0; 
#              all reached minimum library size

#inocas2023Raw <- inocas2023Raw[ !(inocas2023Raw$sample %in% removedSites) , ]

# stopifnot(length(unique(inocas2023Raw$sample[inocas2023Raw$sampleTotalAbd <= minlibrarySize])) == 0)
# 
# print(removedSites)
# remove(removedSites)

# Remove ASVs that don't meet criteria.

# LM 2023-09-26: Added a check for preservation method.


inocas2023Raw <- inocas2023Raw[ inocas2023Raw$primerExpectedLength == "in range", ]
# 2023-09-26 - Result: it dropped 6717 lines; now 32461

# 2023-09-26 - Results: It dropped 29482; now 2979, but all the sites
#              reported arthropods (we have data for th 23 samples)
inocas2023Raw <- inocas2023Raw[ inocas2023Raw$phylumBLASTn %in% phylum, ]

inocasASV <- inocas2023Raw

# What is this for? Where is it used?
# Species abundance (absolute) per sample
inocasMatrix <- ez.matrify(inocasASV,
                           species.name = "ASVHeader",
                           site.name = "sample",
                           abundance = "asvAbsoluteAbundance")

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

remove(inocas2023Raw)

# to do list:
#   + add column names for 2023
#   + import 2023 data (even if wrong) to build out
#   + need to look at which ASVs should be combined or not...






