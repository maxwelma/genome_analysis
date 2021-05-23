directory <- "C:/Users/matil/OneDrive/Dokument/htseq"
#metadata_file <- "C:/Users/matil/Downloads/SraRunTable (1).txt"

sampleFiles = list.files(directory)
sampleCondition = c("Embryonic Forelimbs","Embryonic Forelimbs","Embryonic Hindlimbs","Embryonic Hindlimbs") #condition = c("SRR1719014", "SRR1719015", "SRR1719016", "SRR1719017") c("Embryonic Forelimbs","Embryonic Forelimbs","Embryonic Hindlimbs","Embryonic Hindlimbs")

sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition)

library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~condition)
ddsHTSeq



dds <- DESeq(ddsHTSeq)
res <- results(dds)
res

res_sign_p = subset(res,pvalue<0.05)
res_sign_p@rownames #names of genes with a p-value lower than 0.05

plotMA(res, ylim=c(-7,7))


d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

dds2 = rlog(dds)#Regularized log transformation

#col_means = colMeans(assay(dds2))
#htseq_14 = assay(dds2)[,1]-col_means[1]
#htseq_15 = assay(dds2)[,1]-col_means[2]
#htseq_16 = assay(dds2)[,1]-col_means[3]
#htseq_17 = assay(dds2)[,1]-col_means[4]
#dds_a = cbind(htseq_14,htseq_15,htseq_16,htseq_17) #mean subtraction
#dds2 = rlog(dds)
#pheatmap(assay(dds2),scale = 'row')
#pheatmap(assay(dds2))
#pheatmap(dds_a)


heatmap.2(assay(dds2),scale = 'row',margins=c(17,8))
heatmap.2(assay(dds2[res_sign_p@rownames]),scale = 'row',margins=c(17,8))
par(mar=c(5,6,4,1)+.1)

vsd <- varianceStabilizingTransformation(dds,)
plotPCA(vsd, intgroup=c("condition"))

library("apeglm") #Bayesian shrinkage estimators for effect sizes

resLFC <- lfcShrink(dds, coef="condition_Embryonic.Hindlimbs_vs_Embryonic.Forelimbs", type="apeglm") #Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes.
resLFC_sign <- lfcShrink(dds[res_sign_p@rownames], coef="condition_Embryonic.Hindlimbs_vs_Embryonic.Forelimbs", type="apeglm") #Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes.


plotMA(resLFC, ylim=c(-1,1)) #plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet
plotMA(resLFC_sign, ylim=c(-1,1))





