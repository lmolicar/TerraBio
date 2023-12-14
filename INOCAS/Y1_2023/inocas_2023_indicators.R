#
# inocas_2023_data_indicators
# ***************************
#
# Purpose:
#          This analysis is based on soil samples collected at INOCAS and proc-
#          essed by EcoMol.
#          Takes the eDNA species table provided by EcoMol and determines the
#          TerraBio biodiversity indicators. There are also some species accumu-
#          lation curves to verify sampling completeness.
#
#
# Date             Written by                 Description of changes
# ***********  *****************  **********************************************
# Sep 2023     Karen Dyson        Original code.
#
# Dec 14 2023  Luis Molina        Adapted to read INOCAS data 2023.
#
#


## ----- Data ingestion & setup -----------------------------------

## ---- Importing libraries ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(stringr)
library(compositions)
library(zCompositions)
library(iNEXT)

source("../../allianceBranding.R")

source("../../functions.R")

source("../../../RCode/KDyson_R_Scripts/triplet_fixer.R") # From Karen's github repository

source("../../../RCode/KDyson_R_Scripts/repeat_multipatt.R") # ditto

source("inocas_2023_data_processing.R")


rare = 50

# Per treatment
inocasCompTreatment <- compMatrix(inputMatrix = inocasMatrixLetter)

# Per sample
inocasCompMatrix <- compMatrix(inputMatrix = inocasMatrix)


## ----- Accumulation curves ----------------------------------

# create species accumulation curves
plot(specaccum(inocasMatrix[ grepl(pattern = "-CF-R1|-CF-R2|-CF-R3",
                                   x = rownames(inocasMatrix)), ]),
     xlab = "Number of replicates",
     ylab = "Number of species", ylim = c(10,1200))
plot(specaccum(inocasMatrix[ grepl(pattern = "-F-R1|-F-R2|-F-R3",
                                   x = rownames(inocasMatrix)), ]),
     add = TRUE, col = "green")
plot(specaccum(inocasMatrix[ grepl(pattern = "-I-R1|-I-R2|-I-R3",
                                   x = rownames(inocasMatrix)), ]),
     add = TRUE, col = "blue")
# plot(specaccum(hortaMatrix[ grepl(pattern = "R1-Sy1|R1-Sy2|R1-Sy3",
#                                    x = rownames(hortaMatrix)), ]),
#      add = TRUE, col = "red")
legend(x = 0.5, y = 1300,
       legend = c("Counterfactual", "Control", "Intervention"),
       fill = c("black", "green", "blue"), bty = "n", bg = "n",
       cex = 1)


# create species accumulation curves (remove rare species)
plot(specaccum(inocasMatrix[ grepl(pattern = "-CF-R1|-CF-R2|-CF-R3",
                                   x = rownames(inocasMatrix)), 
                             colSums(inocasMatrix) > rare]),
     xlab = "Number of replicates [rare = 50]",
     ylab = "Number of species", ylim = c(10,350))
plot(specaccum(inocasMatrix[ grepl(pattern = "-F-R1|-F-R2|-F-R3",
                                  x = rownames(inocasMatrix)),
                            colSums(inocasMatrix) > rare]),
     add = TRUE, col = "green")
plot(specaccum(inocasMatrix[ grepl(pattern = "-I-R1|-I-R2|-I-R3",
                                  x = rownames(inocasMatrix)),
                            colSums(inocasMatrix) > rare]),
     add = TRUE, col = "blue")
legend(x = 0.75, y = 375,
       legend = c("Counterfactual", "Control", "Intervention"),
       fill = c("black", "green", "blue"), bty = "n", bg = "n",
       cex = 1)



## ----- Proposed Indicator 1: Alpha diversity  --------------------------------

# Create a table with the alpha diversity measures for each replicate. 
## *** NEED TO ADD FIELD (REPLICATE)***

groupsLUT <- data.frame(lucode = c("CF", "F", "I"),
                        luname = c("Counterfactual",
                                   "Forest",
                                   "Intervention"))
 
inocasGroups <- groupsLUT[match(str_split(rownames(inocasMatrix), "-", simplify = TRUE)[,3], groupsLUT), 2]

groupNames <- str_split(rownames(inocasMatrix), "-", simplify = T)[,1:3]

groupNames <- paste(groupNames[,1], groupNames[,2], groupNames[,3], sep = "-")

inocasAlphaGroup <- alphaGroupMetrics(inocasMatrix, groupNames = groupNames)
 
 inocasAlphaSample <- inocasAlphaGroup %>% mutate(treat = str_split(inocasAlphaGroup$siteType, "-", simplify = T)[,3]) %>% left_join(., y = groupsLUT, by= c("treat" = "lucode"))


# I left Karen's code here to read and understand later ****************** Begin
# Test to compare site diversities between land use types

# Test for raw species number
# 
# speciesRLMER <- lme4::lmer( speciesRichness ~
#                                   siteType +
#                                   (1 | siteType),
#                               data = hortaAlpha,
#                               REML = TRUE)
# anova(speciesRLMER)
# lmerTest::rand(speciesRLMER)
# summary(speciesRLMER)
# test_speciesRLMER<-car::Anova(mod = speciesRLMER)
# emmeans::emmeans(speciesRLMER, pairwise~siteType)    
# ************************************************************************** End



inocasAlphaSample %>%

ggplot(aes(luname, effectiveSR)) + geom_boxplot() + geom_jitter(aes(col=luname, size = I(2)), width = 0.03) +
  labs(color = "Site Type", y = "Effective Species Richness",
         x = "Site Type") +
  scale_color_manual(values = supportingColorPalette) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14))

ggsave("inocas2023_esr.png", width = 8, height = 5, units = "in", dpi = 300)


# Test - After applying a square root transformation
inocasAlphaGroup <- alphaGroupMetrics(sqrt(inocasMatrix),
                                      groupNames = groupNames)



inocasAlphaSample <- inocasAlphaGroup %>% mutate(treat = str_split(inocasAlphaGroup$siteType, "-", simplify = T)[,3]) %>% left_join(., y = groupsLUT, by= c("treat" = "lucode"))


inocasAlphaSample %>%

ggplot(aes(luname, effectiveSR)) + geom_boxplot() + geom_jitter(aes(col=luname, size = I(2)), width = 0.03) +
  labs(color = "Site Type", y = "Effective Species Richness",
       x = "Site Type",
       caption = "After applying square root transform.") +
  scale_color_manual(values = supportingColorPalette) + 
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14))

ggsave("inocas2023_esr_sqrt.png", width = 8, height = 5, units = "in", dpi = 300)


## ----- Proposed Indicator 2: Beta diversity w/ Aitchison distance

# Aitchison distance uses the euclidian distance of the compositional data that
# has been center log transformed; see
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7811025/



    

## ----- PI3: Change in beta diversity -----------------------------------------

# Distance between sites (samples)

inocasASV$sampleID <- paste0("Site",
                             inocasASV$metadata_2, "-",
                             groupsLUT$luname[match(inocasASV$metadata_3, groupsLUT$lucode)])

inocasASV$treatment <- groupsLUT$luname[match(inocasASV$metadata_3,  groupsLUT$lucode)]

inocasSample <- inocasASV %>%
  dplyr::select(sampleID, ASVHeader, asvAbsoluteAbundance) %>%
  group_by(sampleID, ASVHeader) %>%
  summarise(abundance = sum(asvAbsoluteAbundance))

inocasMatrixSample <- ez.matrify(inocasSample, species.name = "ASVHeader", site.name = "sampleID", abundance = "abundance")

# Per sample
inocasCompMatrix <- compMatrix(inputMatrix = inocasMatrixSample)

remap <- inocasASV %>% 
  transmute(original = inocasASV$sampleID,
            pretty = inocasASV$sampleID,
            treatment = inocasASV$treatment) %>% distinct()

levelOrder = c(
  "Counterfactual-Counterfactual",
  "Forest-Forest",
  "Intervention-Intervention",
  "Counterfactual-Forest",
  "Counterfactual-Intervention",
  "Forest-Intervention"
)

aitchisonSample <- vegdist(inocasCompMatrix, "euc", diag = FALSE)

min(aitchisonSample)
max(aitchisonSample)

ait_comparison_sample <- aitComparison(
  inputDist = aitchisonSample,
  remap = remap,
  repeatSamples = TRUE,
  levelsPlot = levelOrder
) + 
  theme(axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14)) +
  labs(title = "Distance between sites") +
  theme(plot.title = element_text(hjust = 0.5))

  ait_comparison_sample
  
  ggsave("ait_inocas_sample_box.png", width = 8, height = 7, units = "in", dpi = 300)
  

# Looking at the overall distance between treatments: heatmap

rownames(inocasCompTreatment) <- groupsLUT$luname[match(rownames(inocasCompTreatment), groupsLUT$lucode)]

aitTreatment <- vegdist(inocasCompTreatment, "euc")

treatmentHeatmap <- aitHeatmap(aitTreatment,
                           fillColor1 = supportingColorPalette[2],
                           fillColor2 = corporateColorPalette[4])
  
treatmentHeatmap + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))

ggsave("inocas_ait_Heatmap_treatment.png", width = 8, height = 7, units = "in", dpi = 300)

## ----- PI 4: qualitative assessment ------------------------------------------
library("FactoMineR")
library("factoextra")
# Good tutorial here: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

# create plot pcas
pca_plots <- inocasCompMatrix %>%
    PCA(., scale.unit = F, graph = F)

viz_pcaPlots <- fviz_pca_ind(
    pca_plots,
    geom.ind = "point",
    col.ind = str_split(rownames(inocasCompMatrix), pattern = "-", simplify = T)[,2],
    addEllipses = T,
    ellipse.type = "confidence",
    pointsize = 3,
    mean.point = F,
    legend.title = "Site Type",
    palette = supportingColorPalette[c(1,3,4,2)]
)
ggpubr::ggpar(viz_pcaPlots,
              title = paste0("Community Composition Visualization using PCA"),
              subtitle = paste0(phylum, collapse = " "), xlab = F, ylab = F, tickslab = F
              )
ggsave("HortaPCA_2022.pdf",
       plot = last_plot(),
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 5,
       units = "in",
       dpi = 300
)


viz_pcaPlots_contrib <- fviz_contrib(pca_plots, choice = "ind", axes = 1:2)

fviz_pca_biplot(pca_plots,
                # Sites
                col.ind = hortaGroups,
                addEllipses = T,
                ellipse.type = "convex",
                label = "var",
                repel = T,
                max.overlaps = 5,
                alpha.var ="contrib")







# I have kept Karen's code here to read and understand later ************* Begin
# create plot pcas for both 
# pca_plotsAll <- hortaCompR1 %>%
#     PCA(., scale.unit = F, graph = F)
# 
# viz_pcaPlotsAll <- fviz_pca_ind(
#     pca_plotsAll,
#     geom.ind = "point",
#     col.ind = c(paste0(hortaGroups, "-b"), paste0(hortaGroups, "-s")),
#     addEllipses = T,
#     ellipse.type = "convex",
#     legend.title = "Group"
# )
# ggpubr::ggpar(viz_pcaPlotsAll,
#               title = paste0("Plots - PCA - ", paste0(phylum, collapse = " "), " - R1"))
# ***************************************************************************End
