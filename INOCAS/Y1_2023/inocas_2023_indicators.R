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
# Dec 19 2023  Karen Dyson        Comments added/code review


## ----- Data ingestion & setup -----------------------------------

## ---- Importing libraries ----
library(iNEXT)

source("../../allianceBranding.R")

source("../../functions.R")

source("../../../RCode/KDyson_R_Scripts/triplet_fixer.R") # From Karen's github repository. Path adjusted to local access.

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

# A lookup table to translate the land use code into a "nice" name.
groupsLUT <- data.frame(lucode = c("CF", "F", "I"),
                        luname = c("Counterfactual",
                                   "Forest",
                                   "Intervention"))

# KD 2023-09-19: This next line of code doesn't work.

# LM 2024-01-11: Fixed.
 
sampleTreatments <- groupsLUT[match(str_split(rownames(inocasMatrix), "-", simplify = TRUE)[,3], groupsLUT$lucode), 2]

# KD 2023-09-19: Why is this recreating the rownames? breaking apart and gluing back together...
# LM 2024-01-11: I wanted to build sample names by dropping the replicate code. Is the purpose clearer this way? I also renamed the output variable
sampleNames <- unlist(lapply(rownames(inocasMatrix),
                             FUN = function(x){str_sub(x, 1, tail(unlist(gregexpr("-", x)), 1)-1)}))

#sampleNames <- str_split(rownames(inocasMatrix), "-", simplify = T)[,1:3]

#groupNames <- paste(groupNames[,1], groupNames[,2], groupNames[,3], sep = "-")

# KD 2023-09-19: No metrics calculated for the replicates? Looks like only the sites, not the site types either. This is the main one used, but the others can help with interpretation and diagnostics.
# LM 2024-01-11: Thanks for advice. I will include diagnostics in a new version after this round of checks.

#inocasSampleAlpha <- alphaGroupMetrics(inocasMatrix, groupNames = groupNames)

# KD 2023-09-19: Not sure what inocasAlphaSample is trying to do--"sample" implies something different than "group" but it's joining two tables.

#inocasAlphaSample <- inocasAlphaGroup %>% mutate(treat = str_split(inocasAlphaGroup$siteType, "-", simplify = T)[,3]) %>% left_join(., y = groupsLUT, by= c("treat" = "lucode"))

# LM 2024-01-11: I was not clear about it. I wanted to estimate ESR for each sample/site and plot them by treatment. Is the purpose clearer now?
 
inocasSampleAlpha <- alphaGroupMetrics(inocasMatrix, groupNames = sampleNames) %>% mutate(treat = str_split(inocasAlphaGroup$siteType, "-", simplify = T)[,3]) %>% left_join(., y = groupsLUT, by= c("treat" = "lucode"))


# KD 2023-09-19: This code tests for significance using LMER. 
 
 
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


# KD 2023-09-19: Note colors don't fit the standards--green should be forest etc. There's a ppt that details this.
# LM 2024-01-11: I think I fixed it. 
 
inocasSampleAlpha %>%

ggplot(aes(luname, effectiveSR)) + geom_boxplot() + geom_jitter(aes(col=luname, size = I(2)), width = 0.03) +
  labs(color = "Site Type", y = "Effective Species Richness",
         x = "Site Type") +
  scale_color_manual(values = supportingColorPalette[c(1,3,4)]) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14))

ggsave("inocas2023_esr.png", width = 8, height = 5, units = "in", dpi = 300)


# KD 2023-09-19: Suggest not overwriting an existing table with new values when they contain different data. It is not super impactful here but as a general rule.

# LM 2024-01-11: Thanks for the advice. Now generating a new table for this.

# Test - After applying a square root transformation
inocasSampleSqrtAlpha <- alphaGroupMetrics(sqrt(inocasMatrix), groupNames = sampleNames) %>% mutate(treat = str_split(inocasAlphaGroup$siteType, "-", simplify = T)[,3]) %>% left_join(., y = groupsLUT, by= c("treat" = "lucode"))



#inocasAlphaSample <- inocasAlphaGroup %>% mutate(treat = str_split(inocasAlphaGroup$siteType, "-", simplify = T)[,3]) %>% left_join(., y = groupsLUT, by= c("treat" = "lucode"))


inocasSampleSqrtAlpha %>%

ggplot(aes(luname, effectiveSR)) + geom_boxplot() + geom_jitter(aes(col=luname, size = I(2)), width = 0.03) +
  labs(color = "Site Type", y = "Effective Species Richness",
       x = "Site Type",
       caption = "After applying square root transform.") +
  scale_color_manual(values = supportingColorPalette[c(1,3,4)]) + 
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


# LM 2024-01-11: Code not required ********************************** Begin
#                it repeats what was done during data processing

inocas2023Filtered$sampleID <- paste0("Site",
                                     inocas2023Filtered$metadata_2, "-",
                            groupsLUT$luname[match(inocas2023Filtered$metadata_3, groupsLUT$lucode)])

inocas2023Filtered$treatment <- groupsLUT$luname[match(inocas2023Filtered$metadata_3,  groupsLUT$lucode)]

# LM 2024-01-12: I haven't done any filtering, so I think I do not need to summarize. By specifying the sampleID the function ez.matrify will aggregate the values per sample. Am I right?

# inocasSample <- inocas2023Filtered %>%
#    dplyr::select(sampleID, ASVHeader, asvAbsoluteAbundance) %>%
#    group_by(sampleID, ASVHeader) %>%
#    summarise(abundance = sum(asvAbsoluteAbundance))

#inocasMatrixSample <- ez.matrify(inocasSample, species.name = "ASVHeader", site.name = "sampleID", abundance = "abundance")

inocasMatrixSample <- ez.matrify(inocas2023Filtered, species.name = "ASVHeader", site.name = "sampleID", abundance = "asvAbsoluteAbundance")

# Per sample
inocasCompMatrix <- compMatrix(inputMatrix = inocasMatrixSample)

# ******************************************************************* End

remap <- inocas2023Filtered %>% 
  transmute(original = inocas2023Filtered$sampleID,
            pretty = inocas2023Filtered$sampleID,
            treatment = inocas2023Filtered$treatment) %>% distinct()

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
    palette = supportingColorPalette[c(1,3,4)]
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


# KD 2023-09-19: Suggest also making the replicate PCAs.

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





# KD 2023-09-19: This code just compares the buffer and silica.

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
