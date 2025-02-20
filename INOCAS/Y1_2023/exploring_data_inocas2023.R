lulut <- data.frame(code = c("CF", "F", "I"),
                    lutype = c("Counterfactual", "Floresta", "Intervention"))


# KD 2023-09-19: Why do these functions have no arguments? Also do not recommend
# this approach to creating tables or doing it like this within a function.
# lulut should be created outside the function and passed as an argument.
# Otherwise you have to change the function for each different dataset...
# in additon, some tables are not passed to this function, so this function will not
# work. The function has no idea what inocas2023Info is.

# qc_plot_purity <- function(){
#   
#   lulut <- data.frame(rbind(
#     c("CF", "Counterfactual"),
#     c("F",  "Floresta"),
#     c("I",  "Intervention")
#   ))
#   
#   colnames(lulut) <- c("code", "lutype")
#   
#   # Purity was imported as character -> convert to numeric
#   inocas2023Info <- inocas2023Info%>% mutate(purityDNA = as.numeric(str_trim(purityDNA)))
#   
#   # KD 2023-09-19: This was imported as character due to white space in the input table. 
#   
#   
#   inocas2023Info %>% select(EcoMolID, purityDNA) %>% 
#     mutate(code = str_split(EcoMolID, "_", simplify = T)[,4]) %>%
#     left_join(., y = lulut, by = "code") %>%
#     ggplot(aes(x = lutype, y = purityDNA)) + geom_boxplot() + 
#     geom_jitter() + 
#     labs(x = "Land use type", y = "Purity (nm)")
# }

library(dplyr)
library(ggplot2)
library(stringr)


# The lookup table should have two columns named code and treatment
Horta2024lut <- data.frame(code = c("Sy", "CF", "F", "R"),
                           treatment = c("Syntropic", "Counterfactual",
                                         "Forest", "Restoration"))
path_sampleInfo <- "h:/My Drive/CAL-PSE Google Drive/04. CAL-PSE Evaluation Activities/TERRABIO/B. TECHNICAL/b. TerraBio_Biodiversity_Component/terrabio_biodiversity_indicators/Horta/2024/HortaSamples2024.csv"


qc_plot_purity <- function(filename, lulut){
  
  df <- read.csv(filename)
  
  # Change the column names
  colnames(df) <- c("Number",	"sampleID",	"treatment",	"labVolume_ml",	"amplification_success",	"replicate",	"EcoMolID",	"concentrationDNA_nguL",	"purityDNA")
  
  # Purity was imported as character -> convert to numeric
  sampleInfo <- df %>% mutate(purityDNA = as.numeric(str_trim(purityDNA)))
  
  # KD 2023-09-19: This was imported as character due to white space in the input table. 
  
  
  sampleInfo %>% select(sampleID, purityDNA) %>% 
    mutate(code = str_split_i(sampleID, "-", 5)) %>%
    left_join(., y = lulut, by = "code") %>%
    ggplot(aes(x = treatment, y = purityDNA)) + geom_boxplot() + 
    geom_jitter() + 
    labs(x = "Treatment", y = "Purity (nm)")
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



inocas2023Info$concentrationDNA_nguL <- as.numeric(str_trim(inocas2023Info$concentrationDNA_nguL))

inocas2023Info$purityDNA <- as.numeric(str_trim(inocas2023Info$purityDNA))

lut_landuse <- data.frame(rbind(c("CF", "Counterfactual"),
                     c("I", "Intervention"),
                     c("F", "Reference")))
colnames(lut_landuse) <- c("lucode","ludesc")


inocas2023Info <- inocas2023Info %>% mutate(lucode = str_split(inocas2023Info$storage, "_", simplify = TRUE)[,3]) %>% left_join(.,lut_landuse, by = "lucode")

# ----- DNA Concentration

ggplot(inocas2023Info, aes(x = ludesc, y = concentrationDNA_nguL)) + geom_boxplot() + ylab( "Concentration ("~ng~mu~L^{-1}~")" ) +
  xlab("Land use") +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 12)) +
  ggtitle("Sample DNA concentration") +
  theme(plot.title = element_text(hjust = 0.5))

outdir <- "d:/OneDrive - CGIAR/Documents/CALPSE/Terrabio/ABF_DEALS/INOCAS/BIODIVERSITY_INDICATORS/ANALYSIS/IMAGES"
filename <- "concentration_per_landuse.png"
dst <- paste(outdir, filename, sep = "/")
ggsave(filename = dst, width = 6, height = 4, units = "in", dpi = 300)


### DNA purity

ggplot(inocas2023Info, aes(x = ludesc, y = purityDNA)) + geom_boxplot() + ylab( expression("A"[260]~"/A"[280])) +
  xlab("Land use") +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 12)) +
  ggtitle("DNA purity") +
  theme(plot.title = element_text(hjust = 0.5))

outdir <- "d:/OneDrive - CGIAR/Documents/CALPSE/Terrabio/ABF_DEALS/INOCAS/BIODIVERSITY_INDICATORS/ANALYSIS/IMAGES"
filename <- "purity_per_landuse.png"
dst <- paste(outdir, filename, sep = "/")
ggsave(filename = dst, width = 6, height = 4, units = "in", dpi = 300)

# ---- ASV count by phylum (unfiltered data) ----
inocas2023Raw %>% select(phylumBLASTn) %>%
  ggplot(aes(x = phylumBLASTn)) + geom_bar() + coord_flip() + 
  xlab("Phylum") + ylab("ASV Count")
outdir <- "d:/OneDrive - CGIAR/Documents/CALPSE/Terrabio/ABF_DEALS/INOCAS/BIODIVERSITY_INDICATORS/ANALYSIS/IMAGES"
filename <- "asv_count_phylum_unfiltered_v1.png"
dst <- paste(outdir, filename, sep = "/")
ggsave(filename = dst, width = 6, height = 6, units = "in", dpi = 300)

# ---- ASV count by phylum (unfiltered data)

dummy <- inocas2023Raw %>% select(phylumBLASTn, order_BLASTn) %>%
  #mutate(code = str_split(sample, "-", simplify = TRUE)[,3]) %>%
  filter(phylumBLASTn == "Arthropoda") %>% select(order_BLASTn) %>% table() %>% as.data.frame() %>% arrange(desc(Freq))

dummy %>% ggplot(aes(reorder(order_BLASTn, -Freq), Freq)) + geom_col() + coord_flip() + 
  labs(x = "Order",
       y = "ASV Count")

ggsave("asv_count_order_unfiltered.png", width = 6, height = 6, units = "in", dpi = 300)



#--- Numerous unindentified ASVs; let's remove them ----
inocasRaw2 <- inocas2023Raw %>% filter(phylumBLASTn != "")
inocasRaw2 %>% ggplot(aes(x = phylumBLASTn)) +
  geom_bar() + coord_flip() + xlab("Phylum") + 
  ylab("ASV Count")
outdir <- "d:/OneDrive - CGIAR/Documents/CALPSE/Terrabio/ABF_DEALS/INOCAS/BIODIVERSITY_INDICATORS/ANALYSIS/IMAGES"
filename <- "asv_count_phylum_unfiltered_v2.png"
dst <- paste(outdir, filename, sep = "/")
ggsave(filename = dst, width = 6, height = 6, units = "in", dpi = 300)


### ASV abundance (filtered data)

# Remove columns we don't need for analysis.
head(inocas2023Raw)
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

# 2023-09-26 - Results: It dropped 29482; now 2979, but all the sites
#              reported arthropods (we have data for th 23 samples)
inocas2023Raw <- inocas2023Raw[ inocas2023Raw$phylumBLASTn %in% phylum, ]

#--- ASV Absolute abundance

inocas2023Raw %>% select(metadata_3, asvAbsoluteAbundance) %>%
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

outdir <- "d:/OneDrive - CGIAR/Documents/CALPSE/Terrabio/ABF_DEALS/INOCAS/BIODIVERSITY_INDICATORS/ANALYSIS/IMAGES"
filename <- "ASV-abs-abundance_per_landuse.png"
dst <- paste(outdir, filename, sep = "/")
ggsave(filename = dst, width = 6, height = 4, units = "in", dpi = 300)

###-- which are the ones above 5000 reads
inocasAbove5k <- inocas2023Raw %>% filter(asvAbsoluteAbundance > 5000)


#--- Sample total abundance

inocas2023Raw %>%
  group_by(sample) %>%
  summarise(abundance = sum(asvAbsoluteAbundance)) %>%
  mutate(lucode = str_split(sample, "-", simplify = T)[,3]) %>%
  left_join(., y = lut_landuse, by = "lucode") %>%
  ggplot(aes(x = ludesc, y = abundance)) + geom_point() +
  labs(x = "Land use type", y = "Sample total abundance",
       caption = "Filtered data") +
  ggtitle("Sample total abundance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Abundance") +
  ylab("Land use") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14))

outdir <- "d:/OneDrive - CGIAR/Documents/CALPSE/Terrabio/ABF_DEALS/INOCAS/BIODIVERSITY_INDICATORS/ANALYSIS/IMAGES"
filename <- "sample-total-abundance_per_landuse.png"
dst <- paste(outdir, filename, sep = "/")
ggsave(filename = dst, width = 6, height = 4, units = "in", dpi = 300)


### Effective species richness



### Checking the impact of the large number of ASVs
#-- Applying square root transformation


