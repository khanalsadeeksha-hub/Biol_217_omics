library('DESeq2')
library('RColorBrewer')
library('gplots')
library('ggplot2')
rawCountTable <- read.table('READemption_analysis_2/output/methanosarcina_gene_quanti_combined/gene_wise_quantifications_combined.csv', skip=1, sep='\t', quote='', comment.char='', colClasses=c(rep('character',10), rep('numeric',4)))
countTable <- round(rawCountTable[,11:length(names(rawCountTable))])
colnames(countTable) <- c('mut_R1','mut_R2','wt_R1','wt_R2')
# Select only the libraries of this species
countTable <- countTable[, c('mut_R1','mut_R2','wt_R1','wt_R2')]
libs <- c('mut_R1','mut_R2','wt_R1','wt_R2')
conds <- c('mut', 'mut', 'wt', 'wt')
reps <- c('1', '2', '1', '2')
samples <- data.frame(row.names=libs, condition=conds, lib=libs, replicate=reps)
dds <- DESeqDataSetFromMatrix(countData=countTable, colData=samples, design=~condition)
dds <- DESeq(dds, betaPrior=TRUE)

# PCA plot
pdf('READemption_analysis_2/output/methanosarcina_deseq/deseq_raw/sample_comparison_pca_heatmap.pdf')
rld <- rlog(dds)
pcaData <- plotPCA(rld, 'condition', intgroup=c('condition', 'replicate'), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, 'percentVar'))
print(ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
geom_point(size=3) +
xlab(paste0('PC1: ',percentVar[1],'% variance')) +
ylab(paste0('PC2: ',percentVar[2],'% variance')) +
coord_fixed())
# Heatmap
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- with(colData(dds), paste(lib, sep=' : '))
hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
heatmap.2(mat, trace='none', col = rev(hmcol), margin=c(13, 13))
comp0 <- results(dds, contrast=c('condition','mut', 'wt'))
write.table(comp0, file='READemption_analysis_2/output/methanosarcina_deseq/deseq_raw/deseq_comp_mut_vs_wt.csv', quote=FALSE, sep='\t')
comp1 <- results(dds, contrast=c('condition','wt', 'mut'))
write.table(comp1, file='READemption_analysis_2/output/methanosarcina_deseq/deseq_raw/deseq_comp_wt_vs_mut.csv', quote=FALSE, sep='\t')
