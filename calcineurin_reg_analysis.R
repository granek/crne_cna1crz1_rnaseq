library("DESeq2")

basedir=Sys.getenv("CNA")
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
resultsNames(ddsHTSeq)[3:6]

# resType <- results(dds, "type_single.read_vs_paired.end")
# head(resType)
##----
# dds <- DESeq(dds)
# res <- results(dds)
# res <- res[order(res$padj),]
# head(res)
# fdrcutoff = 0.05
fdrcutoff = 0.2
KICNA1.res = results(ddsHTSeq,"condition_KI_CNA1_vs_WT")
KICNA1.res = KICNA1.res[order(KICNA1.res$padj),]
table(KICNA1.res$padj < fdrcutoff)
##-------
KICRZ1.res = results(ddsHTSeq,"condition_KI_CRZ1_vs_WT")
KICRZ1.res = KICRZ1.res[order(KICRZ1.res$padj),]
table(KICRZ1.res$padj < fdrcutoff)
##-------
KOcna1.res = results(ddsHTSeq,"condition_KO_cna1_vs_WT")
KOcna1.res = KOcna1.res[order(KOcna1.res$padj),]
table(KOcna1.res$padj < fdrcutoff)
KOcna1Genes = row.names(KOcna1.res[which(KOcna1.res$padj < fdrcutoff),])
##-------
KOcrz1.res = results(ddsHTSeq,"condition_KO_crz1_vs_WT")
KOcrz1.res = KOcrz1.res[order(KOcrz1.res$padj),]
table(KOcrz1.res$padj < fdrcutoff)
KOcrz1Genes = row.names(KOcrz1.res[which(KOcrz1.res$padj < fdrcutoff),])
##-------
length(row.names(KOcna1.res))
length(intersect(KOcna1Genes,KOcrz1Genes))
length(setdiff(KOcna1Genes,KOcrz1Genes))
length(setdiff(KOcrz1Genes,KOcna1Genes))

fileend=paste(fdrcutoff*100,"fdr.csv",sep="")
write.csv(intersect(KOcna1Genes,KOcrz1Genes),file=file.path(outdir,paste("cna1ko_crz1ko_intersect",fileend,sep="_")))
write.csv(setdiff(KOcna1Genes,KOcrz1Genes),file=file.path(outdir,paste("cna1ko_unique",fileend,sep="_")))
write.csv(setdiff(KOcrz1Genes,KOcna1Genes),file=file.path(outdir,paste("crz1ko_unique",fileend,sep="_")))
##-------
get.genes.counts = function(comparison,dds,fdr){
    cur.res = results(dds,comparison)
    # return(row.names(cur.res[which(cur.res$padj < fdr),]))
    return(sum(cur.res$padj < fdr,na.rm = TRUE))
}
mapply(FUN=get.genes.counts,resultsNames(ddsHTSeq),MoreArgs = list(dds=ddsHTSeq,fdr=0.2))
mapply(FUN=get.genes.counts,resultsNames(ddsHTSeq),MoreArgs = list(dds=ddsHTSeq,fdr=0.05))

##------- CRZ1 overexpression ------
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

##------- Normalized Counts Table ------
# head(counts(ddsHTSeq,normalized=TRUE))
# with(colData(ddsHTSeq),paste(row.names(colData(ddsHTSeq)),paste(condition, temp, sep="_"),sep=":"))
ddsHTSeq.counts = counts(ddsHTSeq,normalized=TRUE)
coldat = colData(ddsHTSeq)
colnames(ddsHTSeq.counts) = paste(row.names(coldat),paste(coldat$condition, coldat$temp, sep="_"),sep=":")

fileend=paste(fdrcutoff*100,"fdr_counts.csv",sep="")
write.csv(ddsHTSeq.counts[intersect(KOcna1Genes,KOcrz1Genes),],
          file=file.path(outdir,paste("cna1ko_crz1ko_intersect",fileend,sep="_")))
write.csv(ddsHTSeq.counts[setdiff(KOcna1Genes,KOcrz1Genes),],
          file=file.path(outdir,paste("cna1ko_unique",fileend,sep="_")))
write.csv(ddsHTSeq.counts[setdiff(KOcrz1Genes,KOcna1Genes),],
                          file=file.path(outdir,paste("crz1ko_unique",fileend,sep="_")))

##------- Annotate Gene Lists ------
##------------ Load Annotation --------
annotdir = file.path(Sys.getenv("BREM"),"cneo_hybrid_rnaseq/info")
## wget http://fungalgenomes.org/public/cryptococcus/CryptoDB/product_names/Cneo_H99.AHRD.20131001.tab
ahrd.annot = read.delim(file.path(annotdir,"Cneo_H99.AHRD.20131001.tab"),stringsAsFactors = FALSE)
## http://fungidb.org/fungidb/
fungidb.annot = read.delim(file.path(annotdir,"h99_GenesByTaxon_summary.txt"),
                           stringsAsFactors = FALSE, colClasses = c("character", rep("NULL", 2),"character","NULL"),
                           col.names=c("ID", "organism", "genomeloc", "description","empty"))
annot.df = merge(fungidb.annot, ahrd.annot, by="ID",all = TRUE)
rownames(annot.df) = annot.df$ID
annot.df$ID = NULL
##------------ Extract Genes of Interest --------

fileend=paste(fdrcutoff*100,"fdr_annot.csv",sep="")
write.csv(annot.df[intersect(KOcna1Genes,KOcrz1Genes),],
          file=file.path(outdir,paste("cna1ko_crz1ko_intersect",fileend,sep="_")))
write.csv(annot.df[setdiff(KOcna1Genes,KOcrz1Genes),],
          file=file.path(outdir,paste("cna1ko_unique",fileend,sep="_")))
write.csv(annot.df[setdiff(KOcrz1Genes,KOcna1Genes),],
                          file=file.path(outdir,paste("crz1ko_unique",fileend,sep="_")))
##------- Individual Comparisons ------
# WT_24C <-> KO_crz1_24C
# WT_24C <-> KO_cna1_24C
# KO_cna1_24C <->  KO_crz1_24C

# WT_37C <-> KO_crz1_37C
# WT_37C <-> KO_cna1_37C
# KO_cna1_37C <->  KO_crz1_37C

# KO_cna1_24C <->  KO_cna1_37C
# KO_crz1_24C <->  KO_crz1_37C
