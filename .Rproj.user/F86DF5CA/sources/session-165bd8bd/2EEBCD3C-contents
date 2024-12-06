
# GO analysis 

library(limma)
library(edgeR)
library(Mus.musculus)
organism = 'org.Mm.eg.db'
library(organism, character.only = TRUE)
library(biomaRt)
library(AnnotationDbi)
library(pheatmap)
library(annotables)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(tidyverse)
grcm38 # mouse genes


source("functions/gettop10GO.R")


my_logFC_threshold = 0.2




males <- readRDS("results/Elevel_M_limma_results.RDS")

gettop10GO(males, my_showCategory) %>% 
  mutate(comparison = "Oil - OPFR") -> top10_GOterms_M


write.csv(top10_GOterms_M,"results/results_tables/top10_GOterms_M.csv", row.names = F)





low <- readRDS("results/Elevel_low_limma_results.RDS")

gettop10GO(low, my_showCategory) %>% 
  mutate(comparison = "Oil - OPFR") -> top10_GOterms_low


write.csv(top10_GOterms_low,"results/results_tables/top10_GOterms_low.csv", row.names = F)






high <- readRDS("results/Elevel_high_limma_results.RDS")

gettop10GO(high, my_showCategory) %>% 
  mutate(comparison = "Oil - OPFR") -> top10_GOterms_high


write.csv(top10_GOterms_high,"results/results_tables/top10_GOterms_high.csv", row.names = F)



