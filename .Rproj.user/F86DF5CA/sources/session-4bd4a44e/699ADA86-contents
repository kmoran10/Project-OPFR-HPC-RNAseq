
# DESeq2

library(tidyverse)
library(DESeq2)

cts <- read.csv("rawdata/KT_Counts.csv")
cts <- cts %>% column_to_rownames(.,var = "X")

id <- read.csv("rawdata/kim rnaseq phenodata.csv")
coldata <- id %>% column_to_rownames(.,var = "id")

all(rownames(coldata) == colnames(cts))


### "Quick Start"

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ treatment + sex + treatment:sex)

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name = "treatment_OPFR_vs_Oil")
# or to shrink log fold changes association with condition:
res.s <- lfcShrink(dds, coef="treatment_OPFR_vs_Oil", type="apeglm")
# lots of warning/errors



### DESeq pipeline

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

res <- results(dds, name = "treatment_OPFR_vs_Oil")

res


# Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="treatment_OPFR_vs_Oil", type="apeglm")
resLFC

resOrdered <- res[order(res$pvalue),]
resOrdered


library(IHW)

resIHW <- results(dds, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult



plotMA(res, ylim=c(-2,2))

plotMA(resLFC, ylim=c(-2,2))

dev.off()

#Alternative shrinkage estimators
resNorm <- lfcShrink(dds, coef=2, type="normal")
library(ashr)
resAsh <- lfcShrink(dds, coef=2, type="ashr")


par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")


dev.off()

plotCounts(dds, gene=which.min(res$padj), intgroup="treatment")

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="treatment", 
                returnData=TRUE)
ggplot(d, aes(x=treatment, y=count)) + 
  geom_boxplot() +
  geom_point(position=position_jitter(w=0.1,h=0))


mcols(res)$description


write.csv(as.data.frame(resOrdered), 
          file="results/HPC_bytreatment_DESeq_results.csv")



#### Count data transformations


#variance stabilizing transformations
vsd <- vst(dds, blind=FALSE)
#regularized logarithm
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)


# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))



#### Data quality assessment by sample clustering and visualization

# Heatmap of the count matrix

coldata$treatment <- as.factor(coldata$treatment)

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)
df <- as.data.frame(colData(dds)[,c("treatment","sex")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# not sure what use this is.

sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$group, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# not sure what use this is either..... 




plotPCA(vsd, intgroup=c("treatment"))
plotPCA(vsd, intgroup=c("sex"))
plotPCA(vsd, intgroup=c("E_level"))



mat <- assay(vsd)
mm <- model.matrix(~treatment, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup=c("treatment"))



### Tests of log2 fold change above or below a threshold

par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-2.5,2.5)
resGA <- results(dds, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()




