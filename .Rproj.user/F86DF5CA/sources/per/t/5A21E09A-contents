
#  Oil Comparisons -- cutoff ranges tables


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
library(organism, character.only = TRUE)
library(DOSE)
library(EnhancedVolcano)
library(tidyverse)
grcm38 # mouse genes


my_logFC_threshold = 0.2

y1a <- readRDS("results/Oil_limma_results.RDS")

## M_v_low

y1a$M_v_low %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 243
y1a$M_v_low %>% filter(between(logFC, .2, .75))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 225
y1a$M_v_low %>% filter(between(logFC, .75, 1.5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 14 
y1a$M_v_low %>% filter(between(logFC, 1.5, 3))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 0

y1a$M_v_low %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #  158
y1a$M_v_low %>% filter(between(logFC, -0.75, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 149
y1a$M_v_low %>% filter(between(logFC, -1.5, -0.75)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 7
y1a$M_v_low %>% filter(between(logFC, -3, -1.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 1


#top 10 genes up in M_v_low
y1a$M_v_low %>% filter(P.Value < 0.05) %>% arrange(-logFC) %>% left_join(grcm38, by = "symbol") %>% select(symbol,logFC,P.Value,chr,biotype,description) %>%  head(., 10)
#top 10 genes down in M_v_low
y1a$M_v_low %>% filter(P.Value < 0.05) %>% arrange(logFC) %>% left_join(grcm38, by = "symbol") %>% select(symbol,logFC,P.Value,chr,biotype,description) %>%  head(., 10)





## M_v_high

y1a$M_v_high %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 127
y1a$M_v_high %>% filter(between(logFC, .2, .75))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 122
y1a$M_v_high %>% filter(between(logFC, .75, 1.5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 1 
y1a$M_v_high %>% filter(between(logFC, 1.5, 3))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 0

y1a$M_v_high %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #  118
y1a$M_v_high %>% filter(between(logFC, -0.75, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 117
y1a$M_v_high %>% filter(between(logFC, -1.5, -0.75)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 0
y1a$M_v_high %>% filter(between(logFC, -3, -1.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 0


#top 10 genes up in M_v_high
y1a$M_v_high %>% filter(P.Value < 0.05) %>% arrange(-logFC) %>% left_join(grcm38, by = "symbol") %>% select(symbol,logFC,P.Value,chr,biotype,description) %>%  head(., 10)
#top 10 genes down in M_v_high
y1a$M_v_high %>% filter(P.Value < 0.05) %>% arrange(logFC) %>% left_join(grcm38, by = "symbol") %>% select(symbol,logFC,P.Value,chr,biotype,description) %>%  head(., 10)





## low_v_high

y1a$low_v_high %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 227
y1a$low_v_high %>% filter(between(logFC, .2, .75))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 223
y1a$low_v_high %>% filter(between(logFC, .75, 1.5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 4 
y1a$low_v_high %>% filter(between(logFC, 1.5, 3))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 0

y1a$low_v_high %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #  171
y1a$low_v_high %>% filter(between(logFC, -0.75, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 167
y1a$low_v_high %>% filter(between(logFC, -1.5, -0.75)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 2
y1a$low_v_high %>% filter(between(logFC, -3, -1.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 2


#top 10 genes up in low_v_high
y1a$low_v_high %>% filter(P.Value < 0.05) %>% arrange(-logFC) %>% left_join(grcm38, by = "symbol") %>% select(symbol,logFC,P.Value,chr,biotype,description) %>%  head(., 10)
#top 10 genes down in low_v_high
y1a$low_v_high %>% filter(P.Value < 0.05) %>% arrange(logFC) %>% left_join(grcm38, by = "symbol") %>% select(symbol,logFC,P.Value,chr,biotype,description) %>%  head(., 10)





### Need to make these into pretty table
### need to pull out the top 10 up and down in these into other table








