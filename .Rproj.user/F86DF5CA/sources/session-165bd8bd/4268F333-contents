

# specific genes of interest analysis

####### Analysis of given genes DONE --- none were significantly different. Of those that possibly were (Dgkz, Ptpn5, Grp, and Fkbp5), were not when checked limma eFDR values. 


library(tidyverse)
library(lme4)
library(lmerTest)
library(WGCNA)

## get and clean data ####

data <- read.csv("rawdata/KT_Counts.csv")

phenoData <- read.csv("rawdata/kim rnaseq phenodata.csv")

data[1:10, 1:10]
head(phenoData)

# prepare data

data <- data %>% 
  gather(key = "samples", value = "counts", -X) %>% 
  rename(gene = X) %>% 
  inner_join(., phenoData, by = c("samples" = "id")) %>% 
  select(1, 3, 2) %>% 
  spread(key = "samples", value = "counts") %>% 
  column_to_rownames(var = "gene")


# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

# remove genes that are detected as outliers
data <- data[gsg$goodGenes == TRUE,]

# pca - method 2 for finding outliers 

# exclude outlier samples
samples.to.be.excluded <- c('KT021')
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]
### ***  note in results that we exclude KT021 because of higher variance compared to all other samples 



phenoData <- phenoData %>% 
  column_to_rownames(var = "id")

# exclude outlier samples
colData <- phenoData %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)


all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))
#good



#in r, i have a data frame called "data.subset" and one called "colData" that i want processed with tidyverse when possible. first, current rownames need to be actual columns. for data.subset it needs to be $gene and for colData it needs to be "sample". then I want to remove the rows of "data.subset" when the value in over 75% of the rows across all columns is less that 15. Then, I want to normalize the resulting data frame. Then I want to make the data long format -- currently the columns are all sample names, and I want the rows to be the sample names, with the columns to the names of the current row names. this will make the row names match the row names in "colData." Then join these two data frames in a data frame called gene.data 

# Step 1: Convert rownames to columns
data.subset <- data.subset %>%
  rownames_to_column(var = "gene") # Create a "gene" column from rownames

colData <- colData %>%
  rownames_to_column(var = "sample") # Create a "sample" column from rownames

# Step 2: Remove rows from data.subset where over 75% of values are < 15
data.subset <- data.subset %>%
  filter(rowMeans(select(., -gene) < 15) <= 0.75)

# Step 3: Transform to long format (tidy the data)
data.subset_long <- data.subset %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "value")

# Step 4: Join with colData based on "sample"
gene.data <- data.subset_long %>%
  inner_join(colData, by = "sample")






#### analysis ####

# gotta get all counts, filter bad samples (same as other analyses), filter counts too low, normalize expression, then attach colData to sample names. 
# then can do basically:

# GOI <- df %>% filter(symbol == "GOI") %>% FLIP ORIENTATION
# lmer(gene ~ treatment*E_level + (1|sample)) ---- basically this


## cognition related

# BDNF (Brain-Derived Neurotrophic Factor) -- NONE
gene.data %>% filter(gene == "Bdnf") %>% 
  lm(value ~ treatment*E_level,.) -> bdnflm
summary(bdnflm)
anova(bdnflm)

# CREB (cAMP Response Element-Binding Protein) -- NONE
gene.data %>% filter(gene == "Creb1") %>% 
  lm(value ~ treatment*E_level,.) -> creb1lm
summary(creb1lm)
anova(creb1lm)

# ARC (Activity-Regulated Cytoskeleton-Associated Protein) -- NONE
gene.data %>% filter(gene == "Arc") %>% 
  lm(value ~ treatment*E_level,.) -> arclm
summary(arclm)
anova(arclm)

# NR2B (NMDA Receptor Subunit 2B) -- NONE
gene.data %>% filter(gene == "Grin2b") %>% 
  lm(value ~ treatment*E_level,.) -> grin2blm
summary(grin2blm)
anova(grin2blm)

# Gria1 (AMPA Receptor Subunit 1 / GluR1) -- NONE
gene.data %>% filter(gene == "Gria1") %>% 
  lm(value ~ treatment*E_level,.) -> gria1lm
summary(gria1lm)
anova(gria1lm)

# Syn1 (Synapsin I) -- NONE
gene.data %>% filter(gene == "Syn1") %>% 
  lm(value ~ treatment*E_level,.) -> syn1lm
summary(syn1lm)
anova(syn1lm)

# DGKζ - Dgkz (Diacylglycerol Kinase Zeta) -- Gdkz IS lower in OPFR 
gene.data %>% filter(gene == "Dgkz") %>% 
  lm(value ~ treatment*E_level,.) -> Dgkzlm
summary(Dgkzlm)
anova(Dgkzlm)

gene.data %>% filter(gene == "Dgkz") %>% 
  lm(value ~ treatment+E_level,.) -> Dgkzlm2
summary(Dgkzlm2)


# Shank3 -- NONE
gene.data %>% filter(gene == "Shank3") %>% 
  lm(value ~ treatment*E_level,.) -> Shank3lm
summary(Shank3lm)
anova(Shank3lm)


# CB1R --  Cnr1 (Cannabinoid Receptor 1) -- NONE
gene.data %>% filter(gene == "Cnr1") %>% 
  lm(value ~ treatment*E_level,.) -> Cnr1lm
summary(Cnr1lm)
anova(Cnr1lm)

# Zif268 (Egr1 or Early Growth Response 1) -- NONE
gene.data %>% filter(gene == "Egr1") %>% 
  lm(value ~ treatment*E_level,.) -> Egr1lm
summary(Egr1lm)
anova(Egr1lm)


## stress related 

# CRH -- NONE
gene.data %>% filter(gene == "Crh") %>% 
  lm(value ~ treatment*E_level,.) -> Crh
summary(Crh)
anova(Crh)

# CRHR1 -- NONE
gene.data %>% filter(gene == "Crhr1") %>% 
  lm(value ~ treatment*E_level,.) -> Crhr1
summary(Crhr1)
anova(Crhr1)

gene.data %>% filter(gene == "Crhr1") %>% 
  lm(value ~ treatment+E_level,.) -> Crhr1.2
summary(Crhr1.2)

# CRHR2 -- NONE
gene.data %>% filter(gene == "Crhr2") %>% 
  lm(value ~ treatment*E_level,.) -> Crhr2
summary(Crhr2)
anova(Crhr2)

# PACAP -- Adcyap1 -- NONE
gene.data %>% filter(gene == "Adcyap1") %>% 
  lm(value ~ treatment*E_level,.) -> Adcyap1
summary(Adcyap1)
anova(Adcyap1)

# STEP -- Ptpn5??
gene.data %>% filter(gene == "Ptpn5") %>% 
  lm(value ~ treatment*E_level,.) -> Ptpn5
summary(Ptpn5)
anova(Ptpn5)

gene.data %>% filter(gene == "Ptpn5") %>% 
  lm(value ~ treatment+E_level,.) -> Ptpn5.2
summary(Ptpn5.2)

# UCN -- NO EXPRESSION 
gene.data %>% filter(gene == "Ucn") %>% 
  lm(value ~ treatment*E_level,.) -> Ucn
summary(Ucn)
anova(Ucn)

# AVP -- NONE
gene.data %>% filter(gene == "Avp") %>% 
  lm(value ~ treatment*E_level,.) -> Avp
summary(Avp)
anova(Avp)

# POMC -- NONE
gene.data %>% filter(gene == "Pomc") %>% 
  lm(value ~ treatment*E_level,.) -> Pomc
summary(Pomc)
anova(Pomc)

# NR3C2 -- NONE
gene.data %>% filter(gene == "Nr3c2") %>% 
  lm(value ~ treatment*E_level,.) -> Nr3c2
summary(Nr3c2)
anova(Nr3c2)

# NR3C1 -- NONE
gene.data %>% filter(gene == "Nr3c1") %>% 
  lm(value ~ treatment*E_level,.) -> Nr3c1
summary(Nr3c1)
anova(Nr3c1)

# NPY -- NONE
gene.data %>% filter(gene == "Npy") %>% 
  lm(value ~ treatment*E_level,.) -> Npy
summary(Npy)
anova(Npy)

# GRP
gene.data %>% filter(gene == "Grp") %>% 
  lm(value ~ treatment*E_level,.) -> Grp
summary(Grp)
anova(Grp)

gene.data %>% filter(gene == "Grp") %>% 
  lm(value ~ treatment+E_level,.) -> Grp.2
summary(Grp.2)


# FKBP5
gene.data %>% filter(gene == "Fkbp5") %>% 
  lm(value ~ treatment*E_level,.) -> Fkbp5
summary(Fkbp5)
anova(Fkbp5)

gene.data %>% filter(gene == "Fkbp5") %>% 
  lm(value ~ treatment+E_level,.) -> Fkbp5.2
summary(Fkbp5.2)


#### Compiled potentially significant findings: 
#### (But actually we should just report the relevant limma eFDR when talking about these)
TREAT <- readRDS("results/TREATMENT_limma_results.RDS")
OIL <- readRDS("results/Oil_limma_results.RDS")
OPFR <- readRDS("results/OPFR_limma_results.RDS")

gene.data %>% filter(gene == "Dgkz") %>% 
  lm(value ~ treatment+E_level,.) -> Dgkzlm2
summary(Dgkzlm2) # treatment effect: (b = -4170 ± 1323, n = 29, p < 0.01)

TREAT %>% filter(symbol == "Dgkz") # (eFDR  p = 0.56) -- so, very much not significant? One of these that is happenstance. 



gene.data %>% filter(gene == "Ptpn5") %>% 
  lm(value ~ treatment+E_level,.) -> Ptpn5.2
summary(Ptpn5.2) # treatment effect: (b = -619 ± 314, n = 29, p = 0.060)

TREAT %>% filter(symbol == "Ptpn5") # (eFDR  p = 0.32) -- so, very much not significant? One of these that is happenstance. 



gene.data %>% filter(gene == "Grp") %>% 
  lm(value ~ treatment+E_level,.) -> Grp.2
summary(Grp.2) # sex effect: Males (b = -18.3 ± 7.45, n = 29, p < 0.05)

OIL$M_v_low %>% filter(symbol == "Grp") # (
OIL$M_v_high %>% filter(symbol == "Grp") # (
OPFR$M_v_low %>% filter(symbol == "Grp") # (
OPFR$M_v_high %>% filter(symbol == "Grp") # (
# not in limma results? Must have been excluded for some reason -- since the present p val is rel close to 0.05, though, unlikely that it actually passed the eFDR.




gene.data %>% filter(gene == "Fkbp5") %>% 
  lm(value ~ treatment+E_level,.) -> Fkbp5.2
summary(Fkbp5.2)

TREAT %>% filter(symbol == "Fkbp5") # (eFDR  p = 0.16) -- so, not significant? One of these that is happenstance. 
