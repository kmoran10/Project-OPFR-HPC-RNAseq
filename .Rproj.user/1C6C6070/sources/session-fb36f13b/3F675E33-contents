
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


source("functions/SWAPTEST_gettop10GO.R")


my_logFC_threshold = 0.2

# males %>% filter(symbol == "Cdh11")


males <- readRDS("results/Elevel_M_limma_results.RDS")

SWAPTEST_gettop10GO(males, my_showCategory) %>% 
  mutate(comparison = "Oil - OPFR") -> SWAPTEST_top10_GOterms_M


write.csv(SWAPTEST_top10_GOterms_M,"results/results_tables/SWAPTEST_top10_GOterms_M.csv", row.names = F)





low <- readRDS("results/Elevel_low_limma_results.RDS")

SWAPTEST_gettop10GO(low, my_showCategory) %>% 
  mutate(comparison = "Oil - OPFR") -> SWAPTEST_top10_GOterms_low


write.csv(SWAPTEST_top10_GOterms_low,"results/results_tables/SWAPTEST_top10_GOterms_low.csv", row.names = F)






high <- readRDS("results/Elevel_high_limma_results.RDS")

SWAPTEST_gettop10GO(high, my_showCategory) %>% 
  mutate(comparison = "Oil - OPFR") -> SWAPTEST_top10_GOterms_high


write.csv(SWAPTEST_top10_GOterms_high,"results/results_tables/SWAPTEST_top10_GOterms_high.csv", row.names = F)






###### Top GO-related genes 

### GO Top Genes -- Low E-level
mgenes <- grcm38 # mouse genes

treatgo.le <- read.csv("results/results_tables/SWAPTEST_top10_GOterms_low.csv")

# First, process the treatgo data as before
gene_counts.le <- treatgo.le %>%
  separate_rows(geneID, sep = "/") %>%
  count(direction, geneID, name = "gene_count") %>%
  arrange(direction, desc(gene_count)) %>%
  group_by(direction) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Join with mgenes to get the gene descriptions
result_with_desc.le <- gene_counts.le %>%
  # Join with mgenes using geneID = symbol
  left_join(
    mgenes %>% select(symbol, gene_description = description),
    by = c("geneID" = "symbol")
  ) %>%
  # Now join with original data to get the treatgo descriptions
  left_join(
    treatgo.le %>%
      separate_rows(geneID, sep = "/") %>%
      select(direction, geneID, treatgo_description = Description),
    by = c("direction", "geneID")
  ) %>%
  # Group and combine the treatgo descriptions
  group_by(direction, geneID, gene_count, rank, gene_description) %>%
  summarise(
    treatgo_descriptions = paste(unique(treatgo_description), collapse = "; "),
    .groups = "drop"
  ) %>%
  # Keep only top genes
  filter(rank <= 10) %>%
  arrange(direction, desc(gene_count))

result_with_desc.le


GO_top_genes.le <- result_with_desc.le %>%
  mutate(
    treat_data = map(geneID, ~ low %>% filter(symbol == .x))
  )%>%
  unnest(treat_data, keep_empty = TRUE)

SWAPTEST_GO_top_genes_LE_table <- GO_top_genes.le %>% 
  select(1,2,10,3,5,8,9)

write.csv(SWAPTEST_GO_top_genes_LE_table,"results/results_tables/SWAPTEST_GO_top_genes_LE_table.csv", row.names = F)



# high e-level
treatgo.he <- read.csv("results/results_tables/SWAPTEST_top10_GOterms_high.csv")

# First, process the treatgo data as before
gene_counts.he <- treatgo.he %>%
  separate_rows(geneID, sep = "/") %>%
  count(direction, geneID, name = "gene_count") %>%
  arrange(direction, desc(gene_count)) %>%
  group_by(direction) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Join with mgenes to get the gene descriptions
result_with_desc.he <- gene_counts.he %>%
  # Join with mgenes using geneID = symbol
  left_join(
    mgenes %>% select(symbol, gene_description = description),
    by = c("geneID" = "symbol")
  ) %>%
  # Now join with original data to get the treatgo descriptions
  left_join(
    treatgo.he %>%
      separate_rows(geneID, sep = "/") %>%
      select(direction, geneID, treatgo_description = Description),
    by = c("direction", "geneID")
  ) %>%
  # Group and combine the treatgo descriptions
  group_by(direction, geneID, gene_count, rank, gene_description) %>%
  summarise(
    treatgo_descriptions = paste(unique(treatgo_description), collapse = "; "),
    .groups = "drop"
  ) %>%
  # Keep only top genes
  filter(rank <= 10) %>%
  arrange(direction, desc(gene_count))

result_with_desc.he


GO_top_genes.he <- result_with_desc.he %>%
  mutate(
    treat_data = map(geneID, ~ high %>% filter(symbol == .x))
  )%>%
  unnest(treat_data, keep_empty = TRUE)

SWAPTEST_GO_top_genes_HE_table <- GO_top_genes.he %>% 
  select(1,2,10,3,5,8,9)

write.csv(SWAPTEST_GO_top_genes_HE_table,"results/results_tables/SWAPTEST_GO_top_genes_HE_table.csv", row.names = F)



# Males

treatgo.m <- read.csv("results/results_tables/SWAPTEST_top10_GOterms_M.csv")

# First, process the treatgo data as before
gene_counts.m <- treatgo.m %>%
  separate_rows(geneID, sep = "/") %>%
  count(direction, geneID, name = "gene_count") %>%
  arrange(direction, desc(gene_count)) %>%
  group_by(direction) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Join with mgenes to get the gene descriptions
result_with_desc.m <- gene_counts.m %>%
  # Join with mgenes using geneID = symbol
  left_join(
    mgenes %>% select(symbol, gene_description = description),
    by = c("geneID" = "symbol")
  ) %>%
  # Now join with original data to get the treatgo descriptions
  left_join(
    treatgo.m %>%
      separate_rows(geneID, sep = "/") %>%
      select(direction, geneID, treatgo_description = Description),
    by = c("direction", "geneID")
  ) %>%
  # Group and combine the treatgo descriptions
  group_by(direction, geneID, gene_count, rank, gene_description) %>%
  summarise(
    treatgo_descriptions = paste(unique(treatgo_description), collapse = "; "),
    .groups = "drop"
  ) %>%
  # Keep only top genes
  filter(rank <= 10) %>%
  arrange(direction, desc(gene_count))

result_with_desc.m


GO_top_genes.m <- result_with_desc.m %>%
  mutate(
    treat_data = map(geneID, ~ males %>% filter(symbol == .x))
  )%>%
  unnest(treat_data, keep_empty = TRUE)

SWAPTEST_GO_top_genes_M_table <- GO_top_genes.m %>% 
  select(1,2,10,3,5,8,9)


write.csv(SWAPTEST_GO_top_genes_M_table,"results/results_tables/SWAPTEST_GO_top_genes_M_table.csv", row.names = F)


