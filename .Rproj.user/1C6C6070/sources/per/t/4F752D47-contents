
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

coldata$group <- paste0(coldata$treatment, "_", coldata$E_level)

# removing KT021 as outlier
dlNorm$KT021 <- NULL
coldata <- coldata[-21,]
## REMOVING KT021  since marked as outlier.

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

dge.dl_M <- dge.dl[, dge.dl$samples$group %in% c("Oil_M", "OPFR_M")]
dge.dl_M$samples$group <- droplevels(dge.dl_M$samples$group)
dge.dl_M$samples$group
dge.dl<- dge.dl_M
dge.dl$samples$group


coldata %>%
  filter(E_level == "low") %>% 
  dplyr::select(subject, group) -> var_info

#shortening used data to subset

dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check




coldata$group %>% 
  factor(.,levels = c("Oil_M","OPFR_M")) -> group.dl

design.dl <- model.matrix(~ 0 + group.dl)
colnames(design.dl) -> mycolnames


v.dl = voom(dge.dl, design.dl, plot = F)
vfit.dl = lmFit(v.dl, design.dl)


contr.matrix <- makeContrasts(group.dlOil_M-group.dlOPFR_M,
                              levels = design.dl)


vfit.dl2 <- contrasts.fit(vfit.dl, contr.matrix)

efit.dl2 = eBayes(vfit.dl2)

p.dl.limma2 = efit.dl2[["p.value"]]
head(p.dl.limma2)

saveRDS(v.dl, "results/limma_vdl_Elevel_M.RDS")


### random sampling - don't redo unless needed ####
####### THIS WILL LIKELY NOT WORK SINCE ONLY HAVE FEW SAMPLES.


R = 2500
set.seed(312)

#to store pvalues in
p.dl.rand = vector('list',length = R)

# to store "t" values (coefficients)
p.dl.rand.t = vector('list',length = R)

for( g in 1 : R){
  print(paste("Starting on Permutation", g))
  
  # Randomize the traits
  
  group.dl.rand = sample(group.dl)
  
  # Model
  design.dl.rand = model.matrix(~0 + group.dl.rand)
  colnames(design.dl.rand) <- mycolnames
  
  # Calculate p-values based on randomized traits
  v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
  vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
  
  vfit.dl.rand2 <- contrasts.fit(vfit.dl.rand, contr.matrix)
  
  efit.dl.rand2 = eBayes(vfit.dl.rand2)
  
  p.dl.rand[[g]] = efit.dl.rand2[["p.value"]]
  p.dl.rand.t[[g]] = efit.dl.rand2[["t"]]
  
}


q.dl <- Reduce(`+`, lapply(p.dl.rand, \(x) {
  (x < p.dl.limma2)
}))

head(q.dl)

q.dl = q.dl / R

q.dl = as.data.frame(q.dl)


efit.dl2[["p.value"]] <- q.dl

sum(duplicated(row.names(efit.dl2$coefficients)))



saveRDS(q.dl,("results/limma_vdl_cutoff5_2000_tworand_Elevel_M.RDS"))

#### done random sampling - don't redo unless needed ####


q.dl <- readRDS("results/limma_vdl_cutoff5_2000_tworand_Elevel_M.RDS")


##### Analysis pulling genes out for each contrast 
tmp1 <- contrasts.fit(efit.dl2, coef = 1) # 


#keep in mind these P.Values ARE the eFDRs

Elevel_M_limma_results <- topTable(tmp1, sort.by = "P", n = Inf) %>% 
  rownames_to_column('symbol') %>% 
  dplyr::select(symbol,logFC,P.Value)

Elevel_M_limma_results1 <- Elevel_M_limma_results %>% 
  left_join(grcm38, by = "symbol") %>% 
  filter(!is.na(symbol)) %>% 
  filter(!is.na(entrez)) %>%
  select(symbol,logFC,P.Value,chr,entrez,start,end,biotype,description)



saveRDS(Elevel_M_limma_results1,"results/Elevel_M_limma_results.RDS")


Elevel_M_limma_results1 <- readRDS("results/Elevel_M_limma_results.RDS")

Elevel_M_limma_results1 %>% 
  filter(., P.Value<0.05) %>% 
  filter(., P.Value != 0) %>%
  summarise(.,Up = sum(logFC>0.2),
            Down = sum(logFC<0.2)) %>% 
  mutate(.,Total = Up + Down) 

hist(Elevel_M_limma_results1$logFC)







