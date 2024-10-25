
# limma

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

all(rownames(coldata) == colnames(dlNorm))

###

d = apply(dlNorm, 2, as.numeric)
dim(d)

d0 = DGEList(d, group = coldata$group)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)

dge.dl$samples$group

coldata %>% 
  dplyr::select(subject, treatment) -> var_info

coldata$treatment %>% 
  factor(.,levels = c("Oil","OPFR")) -> group.dl











