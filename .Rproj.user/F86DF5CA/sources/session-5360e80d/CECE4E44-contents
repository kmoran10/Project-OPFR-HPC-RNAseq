###### NEED TO DO ALL THIS FOR TREATMENT - JUST COPIED OVER TEMPLATE FOR NOW.

# pairwise comparisons of normalized expression of important genes 
# selected by large logFC, MEcolor strength, and theoretical importance 

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(CorLevelPlot)
library(gridExtra)
library(clusterProfiler)
library(enrichplot)
library(biomaRt)
library(AnnotationDbi)
library(annotables)
grcm38 <- grcm38
library(ggpubr)
library(tidyverse)


#### making normalized counts - shouldn't need to be done again ####
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


table(gsg$goodGenes)
table(gsg$goodSamples)


# remove genes that are detected as outliers
data <- data[gsg$goodGenes == TRUE,]


# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)
# KT021 definitely outlier.  and maybe KT004


# pca - method 2 for finding outliers 

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
# less clear here -- going to go back to my voom transformed data and look at those PCAs -- if KT021 is the oilier over there too, will remove from this analysis and prior analyses too (damnit)
# ok remove KT021 lol

# exclude outlier samples
samples.to.be.excluded <- c('KT021')
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]
### ***  note in results that we exclude KT021 because of higher variance compared to all other samples 



# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset


phenoData <- phenoData %>% 
  column_to_rownames(var = "id")

# exclude outlier samples
colData <- phenoData %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)

# # fixing column names in colData ---- unneccessary here, all clean already
# names(colData)
# 
# # selecting relevant info
# colData <- colData %>% 
#   select(3, 10, 11, 13, 14, 15, 20, 22)


# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))



# create dds
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1) # not specifying model "because we need this DSeq data set to perform variance stabilizing transformation"


## remove all genes with counts < 15 in more than 75% of samples (29*0.75=22)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 22,]
nrow(dds75) # 14068 genes



# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()




norm.counts2 <- norm.counts
norm.counts2 <- as.data.frame(norm.counts2)
norm.counts3 <- tibble::rownames_to_column(norm.counts2, "sample")
colData2 <- tibble::rownames_to_column(colData, "sample")
normcounts.coldata.TREAT <- cbind(colData2, norm.counts3)

write.csv(normcounts.coldata.TREAT, "results/normcounts_coldata_TREATMENT.csv")


#### shouldn't need to do above again ####

normcounts.coldata.TREATMENT <- read.csv("results/normcounts_coldata_TREATMENT.csv")
normcounts.coldata.TREATMENT$X <- NULL
TREATMENT_limma_results1 <- readRDS("results/TREATMENT_limma_results.RDS")
TREATMENT_MEsalmon_limma <- read.csv("results/results_tables/TREATMENT_MEsalmon_limma.csv") #HIGHER expression in OPFR
TREATMENT_MEsalmon_limma$X <- NULL



### top 50 biggest logFC genes 
TREATMENT_limma_results1 %>% 
  arrange(logFC) %>% 
  head(25)

TREATMENT_limma_results1 %>% 
  arrange(-logFC) %>% 
  head(25)

#visualize some favorites 

### TEMPLATE
# Gabra2 - save as w400 h600 -- remember to replace stat compare means with the limma eFDR. 
ggplot(normcounts.coldata.TREATMENT, aes(treatment, Gabra2, fill=treatment)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 1.5) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("skyblue", "hotpink")) +
  labs(title="X. Gabra2 Gene Expression",x="Treatment", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 


