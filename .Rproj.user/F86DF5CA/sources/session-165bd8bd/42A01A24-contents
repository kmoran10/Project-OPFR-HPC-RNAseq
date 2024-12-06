
# initial cleaning and checks

library(limma)
library(tidyverse)
library(scales)


library(tidyverse)
library(DESeq2)
library(pheatmap)
library(annotables)
library(ggpubr)

source("functions/geom_boxjitter.R")


rawcounts <- read.csv("rawdata/KT_Counts.csv")
bcx1 <- rawcounts %>% column_to_rownames(.,var = "X")

id <- read.csv("rawdata/kim rnaseq phenodata.csv")
idxx <- id %>% column_to_rownames(.,var = "id")

all(rownames(idxx) == colnames(bcx1)) #check they match

#### count checks

colSums(bcx1[,]) %>% 
  as_tibble_row() %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'id') %>% 
  dplyr::rename(genecounts = V1) -> genecounts

genecounts <- genecounts %>%
  full_join(id) 


#########


genecounts %>% 
  ggplot(aes(genecounts)) +
  geom_histogram(bins = 54,alpha =0.5,color = 'grey') +
  theme_classic() 

counts <- genecounts %>%
  ggplot(aes(id,genecounts))+
  geom_bar(stat = 'identity')+
  coord_flip()

ggsave("imgs/braincounts_byid.png",counts, height = 5,width= 5, dpi = 150)

# all quite tight -- well processed.


####

counts2 <- genecounts %>%
  ggplot(aes(treatment,genecounts, color = treatment))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 1,
                 width = 0.5,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.6)) +
  theme_classic()+
  theme(legend.position = "none")
counts2

ggsave("imgs/braincounts_bytreatment.png",counts2, height = 5,width= 5, dpi = 150)


## WOW -- simply less expression overall in OPFR treated animals????? -- double check in on this later. 


counts3 <- genecounts %>%
  ggplot(aes(treatment,genecounts, color = sex))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 1,
                 width = 0.5,
                 jitter.height = 0.02, jitter.width = 0.02, errorbar.draw = TRUE,
                 position = position_dodge(0.8)) +
  theme_classic()
counts3


ggsave("imgs/braincounts_bygroupandROI_ham.png",counts3, height = 5,width= 5, dpi = 150)

### OVERALL less expression in OPFR, but more dramatically impacted in Fs than Ms? -- check interaction here. 




#### counts~treatment

genecounts %>%
  ggplot(aes(treatment,genecounts, fill = treatment))+
  stat_compare_means(method = "t.test", size = 6, label.x = 1.5) +
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 1,
                 width = 0.5,
                 jitter.height = 0.02, jitter.width = 0.02, errorbar.draw = TRUE,
                 position = position_dodge(0.8)) +
  scale_fill_manual(values=c("skyblue", "hotpink")) +
  labs(title="Counts~Treatment",x="Group", y = "Gene Counts") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 


#### counts~treatment*sex

genecounts %>%
  ggplot(aes(treatment,genecounts, fill = sex))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 1,
                 width = 0.5,
                 jitter.height = 0.02, jitter.width = 0.02, errorbar.draw = TRUE,
                 position = position_dodge(0.8)) +
  scale_fill_manual(values=c("hotpink", "skyblue")) +
  labs(title="Counts~Treatment*Sex",x="Group", y = "Gene Counts") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
#        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 


lm1 <- lm(genecounts~treatment*sex, data = genecounts)
summary(lm1)
#no interaction between treatment and sex on raw gene count




# what about E-level?
genecounts %>%
  ggplot(aes(treatment,genecounts, fill = E_level))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 1,
                 width = 0.5,
                 jitter.height = 0.02, jitter.width = 0.02, errorbar.draw = TRUE,
                 position = position_dodge(0.8)) +
  scale_fill_manual(values=c("hotpink", "pink", "skyblue")) +
  labs(title="Counts~Treatment*E_Level",x="Group", y = "Gene Counts") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 

lm2 <- lm(genecounts~treatment*E_level, data = genecounts)
summary(lm2)

lm3 <- lm(genecounts~treatment+E_level, data = genecounts)
summary(lm3)

# while it may look it somewhat, no effect or interaction of sex or e level on raw count