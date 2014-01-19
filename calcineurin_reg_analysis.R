library("DESeq2")
## library("pasilla")

basedir=("/Users/josh/Documents/BioinfCollabs/HeitmanLab/calcineurin_reg")
countdir=file.path(basedir,"counts")
counttab.file=file.path(basedir,"info","calcineurin_sample_table.csv")
outdir=file.path(basedir,"results")
sampleComparePlots = file.path(outdir,"sample_compare_plots.pdf")

sampleData = read.csv(counttab.file,colClasses=c("character","numeric","character","factor","factor"))
sampleData$genotype = factor(sampleData$genotype, levels=c("WT", "KI_CNA1", "KI_CRZ1", "KO_cna1", "KO_crz1"))
sampleTable = transform(sampleData, sampleName=sample_name,fileName=sample_file,condition=genotype)
rownames(sampleTable) = sampleTable$sample_num
sampleTable <- subset(sampleTable, select = c(sampleName,fileName,condition,temp) )

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                           directory = countdir,
                                           design= ~ condition)
#colData(x)$condition <- factor(colData(x)$condition, levels=c("Control","A","B"))
##    directory <- system.file("extdata", package="pasilla", mustWork=TRUE)
##    sampleFiles <- grep("treated",list.files(directory),value=TRUE)
##    sampleCondition <- sub("(.*treated).*","\\1",sampleFiles)
##    sampleTable <- data.frame(sampleName = sampleFiles,
##                              fileName = sampleFiles,
##                              condition = sampleCondition)
##    ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
##                                           directory = directory,
##                                           design= ~ condition)
##    colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
##                                          levels=c("untreated","treated"))
## ddsHTSeq

## design(dds) <- formula(~ type + condition)
design(ddsHTSeq) <- formula(~ temp + condition)
ddsHTSeq <- DESeq(ddsHTSeq)
res <- results(ddsHTSeq)
head(res)
resultsNames(ddsHTSeq)
# resType <- results(dds, "type_single.read_vs_paired.end")
# head(resType)
##----
# dds <- DESeq(dds)
# res <- results(dds)
# res <- res[order(res$padj),]
# head(res)
KICNA1.res = results(ddsHTSeq,"condition_KI_CNA1_vs_WT")
KICNA1.res = KICNA1.res[order(KICNA1.res$padj),]
table(KICNA1.res$padj < 0.2)
##-------
KICRZ1.res = results(ddsHTSeq,"condition_KI_CRZ1_vs_WT")
KICRZ1.res = KICRZ1.res[order(KICRZ1.res$padj),]
table(KICRZ1.res$padj < 0.2)
##------- CRZ1 overexpression
cna1.id="CNAG_04796"
crz1.id="CNAG_00156"
crz1mat = counts(ddsHTSeq,normalized=TRUE)[c(cna1.id,crz1.id),]
colnames(crz1mat) = with(colData(ddsHTSeq),paste(condition, temp, sep=" : "))
pdf(file.path(outdir,"crz1_overexpress.pdf"))
heatmap.2(t(crz1mat), trace="none", col = rev(hmcol), margin=c(7, 7),cexCol=1)
dev.off()
##----- Heatmap of the sample-to-sample distances------
library("RColorBrewer")
library("gplots")
# select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30] 
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

rld <- rlogTransformation(ddsHTSeq, blind=TRUE)
distsRL <- dist(t(assay(rld)))
## A heatmap of this distance matrix gives us an overview over similarities and dissimilarities between samples (Figure 8):
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(ddsHTSeq),
                                       paste(condition, temp, sep=" : "))
pdf(sampleComparePlots)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(7, 7))
plotPCA(rld, intgroup=c("condition", "temp"))
dev.off()
