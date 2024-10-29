
# limma all 6 conditions pairwise

library(limma)
library(Glimma)
library(edgeR)

library(biomaRt)
library(AnnotationDbi)
library(annotables)
grcm38 <- grcm38
library(tidyverse)

###

cts <- read.csv("rawdata/KT_Counts.csv")
dlNorm <- cts %>% column_to_rownames(.,var = "X")
dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]

id <- read.csv("rawdata/kim rnaseq phenodata.csv")
coldata <- id %>% column_to_rownames(.,var = "id")

coldata$group <- paste0(coldata$treatment, "_", coldata$E_level)


all(rownames(coldata) == colnames(dlNorm))

###
