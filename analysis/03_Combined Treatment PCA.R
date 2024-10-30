
# DESeq2 PCA combined group
# to determine variance attributable to different groupp factors
  # likely, treatment will have effect
  # visually test if E_level (and therefore sex) has effect

library(tidyverse)
library(DESeq2)

cts <- read.csv("rawdata/KT_Counts.csv")
cts <- cts %>% column_to_rownames(.,var = "X")

id <- read.csv("rawdata/kim rnaseq phenodata.csv")
coldata1 <- id %>% column_to_rownames(.,var = "id")

coldata1$group <- paste0(coldata1$treatment, "_", coldata1$E_level)
coldata <- coldata1


all(rownames(coldata) == colnames(cts))

### DESeq

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ group)


dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

#removing 0s
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]



vsd <- vst(dds, blind=FALSE)

library("vsn")

meanSdPlot(assay(vsd))

library("RColorBrewer")

plotPCA(vsd, intgroup=c("group"))
plotPCA(vsd, intgroup=c("treatment"))
plotPCA(vsd, intgroup=c("E_level"))


pcaData <- plotPCA(vsd, intgroup=c("group","E_level","treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=group.1, label = rownames(pcaData))) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + geom_text() + labs(title="Combined Treatment")


ggplot(pcaData, aes(PC1, PC2, color=E_level, label = treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + geom_text() + labs(title="Combined Treatment")


ggplot(pcaData, aes(PC1, PC2, color=E_level, label = name)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + geom_text() + labs(title="Combined Treatment")

#KT021 is def an outlier 
#KT004 may be, but retain for now.

#### THE NEXT TIME YOU DO THIS IN YOUR LIFE, REMEMBER TO FILTER OUT OUTLIERS. INCLUDE THE OUTLIER CHECK STEPS FROM WGCNA.

# 
# library(WGCNA)
# library(DESeq2)
# library(GEOquery)
# library(CorLevelPlot)
# library(gridExtra)
# library(clusterProfiler)
# library(enrichplot)
# library(biomaRt)
# library(AnnotationDbi)
# library(annotables)
# grcm38 <- grcm38
# library(tidyverse)
# source("functions/gettop10GO.R")
# 
# 
# allowWGCNAThreads()          # allow multi-threading (optional)
# 
# 
# 
# data <- read.csv("rawdata/KT_Counts.csv")
# 
# phenoData <- read.csv("rawdata/kim rnaseq phenodata.csv")
# 
# data[1:10, 1:10]
# head(phenoData)
# 
# # prepare data
# 
# data <- data %>% 
#   gather(key = "samples", value = "counts", -X) %>% 
#   rename(gene = X) %>% 
#   inner_join(., phenoData, by = c("samples" = "id")) %>% 
#   select(1, 3, 2) %>% 
#   spread(key = "samples", value = "counts") %>% 
#   column_to_rownames(var = "gene")
# 
# 
# # 2. QC - outlier detection ------------------------------------------------
# # detect outlier genes
# 
# gsg <- goodSamplesGenes(t(data))
# summary(gsg)
# gsg$allOK
# 
# 
# table(gsg$goodGenes)
# table(gsg$goodSamples)
# 
# 
# # remove genes that are detected as outliers
# data <- data[gsg$goodGenes == TRUE,]
# 
# 
# # detect outlier samples - hierarchical clustering - method 1
# htree <- hclust(dist(t(data)), method = "average")
# plot(htree)