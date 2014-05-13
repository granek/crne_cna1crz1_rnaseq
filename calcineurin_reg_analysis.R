library("DESeq2")
library("RColorBrewer")
library("gplots")

args <- commandArgs(trailingOnly = TRUE)
print(args)


if ("--usecwd" %in% args) {
    ## outdir<<-file.path(".","results")
    basedir<<-"."
} else {
    ## outdir<<-file.path(basedir,"results")
    basedir<<-Sys.getenv("CNA")
}
outdir=file.path(basedir,"results")
annotdir = file.path(basedir,"info")
countdir=file.path(basedir,"counts")
dir.create(outdir)

counttab.file=file.path(basedir,"info","calcineurin_sample_table.csv")
# sampleComparePlots = file.path(outdir,"sample_compare_plots.pdf")
sampleCompareHeatmap = file.path(outdir,"sample_compare_heatmap.pdf")
sampleComparePCA = file.path(outdir,"sample_compare_pca.pdf")

sampleData = read.csv(counttab.file,colClasses=c("character","numeric","character","factor","factor"))
sampleData$genotype = factor(sampleData$genotype, levels=c("WT", "KI_CNA1", "KI_CRZ1", "KO_cna1", "KO_crz1"))
sampleTable = transform(sampleData, sampleName=sample_name,fileName=sample_file,condition=genotype)
rownames(sampleTable) = sampleTable$sample_num
sampleTable <- subset(sampleTable, select = c(sampleName,fileName,condition,temp) )

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                           directory = countdir,
                                           design= ~ condition)
## design(dds) <- formula(~ type + condition)
design(ddsHTSeq) <- formula(~ temp + condition)
ddsHTSeq <- DESeq(ddsHTSeq)
res <- results(ddsHTSeq)

##===========================================================================
##===========================================================================
## HERE >>>>>>>
## stop("Set up filtering")
## if ("--nofilter" %in% args) {
##     print("Not Filtering")
## } else {
##     print("Filtering Genes with mean counts less than 10 or NA pvalue")
##     use <- res$baseMean >= 10 & !is.na(res$pvalue)
##     table(use)
##     resFilt <- res[use,]
##     resFilt$padj <- p.adjust(resFilt$pvalue, method="BH")
##     sum(res$padj < .1, na.rm=TRUE)
##     sum(resFilt$padj < .1, na.rm=TRUE)
##     res = resFilt
##     ## antisense.res = results(anti.ddsHTSeq)
##     ## antisense.use <- antisense.res$baseMean >= 10
##     ## antisense.counts.filt = antisense.counts[antisense.use,]
## }
## 
## HERE <<<<<<
##===========================================================================
##===========================================================================


head(res)
resultsNames(ddsHTSeq)
resultsNames(ddsHTSeq)[3:6]



FindDiffGenes = function(ddsHTSeq,outdir, fdrcutoff=0.05, fccutoff=2,countfilter=FALSE){
    ## stop("need to do filtering for each comparison????")
    log2fc = log2(fccutoff)
    # for(var in seq) expr
    ListOfGeneVecs = list()
    for(sample in c("KI_CNA1","KI_CRZ1","KO_cna1","KO_crz1")) {
        coeff = paste("condition",sample, "vs_WT",sep="_")
        ## print(coeff)
        cur.res = results(ddsHTSeq,coeff)
        cur.res = cur.res[order(cur.res$padj),]
        ## print(sample)
        print(
            table(cur.res$padj < fdrcutoff,
                  abs(cur.res$log2FoldChange) >= log2fc,
                  dnn=c(paste("FDR<",fdrcutoff), paste("FC>",fccutoff)))
            )
        if (countfilter) {
            print("Filtering Genes with mean counts less than 10 or NA pvalue")
            use <- cur.res$baseMean >= 10 & !is.na(cur.res$pvalue)
            table(use)
            resFilt <- cur.res[use,]
            resFilt$padj <- p.adjust(resFilt$pvalue, method="BH")
            sum(cur.res$padj < .1, na.rm=TRUE)
            sum(resFilt$padj < .1, na.rm=TRUE)
            cur.res = resFilt
            print(
                table(cur.res$padj < fdrcutoff,
                      abs(cur.res$log2FoldChange) >= log2fc,
                      dnn=c(paste("FDR<",fdrcutoff), paste("FC>",fccutoff)))
            )

            filtered="_filt"
        } else {
            filtered=""
        }
        fileend=paste(fdrcutoff*100,"fdr_", fccutoff,"fc",filtered,".csv",sep="")
        filt.res = cur.res[which((cur.res$padj < fdrcutoff) &
                          (abs(cur.res$log2FoldChange) >= log2fc)),]
        write.csv(as.data.frame(filt.res),file=file.path(outdir,paste(coeff, "results",fileend,sep="_")))
        ListOfGeneVecs[[sample]] = row.names(filt.res)
    }

    KOcna1Genes = ListOfGeneVecs[["KO_cna1"]]
    KOcrz1Genes = ListOfGeneVecs[["KO_crz1"]]
    write.csv(KOcna1Genes,file=file.path(outdir,paste("cna1ko_genes",fileend,sep="_")))
    write.csv(KOcrz1Genes,file=file.path(outdir,paste("crz1ko_genes",fileend,sep="_")))
    write.csv(intersect(KOcna1Genes,KOcrz1Genes),file=file.path(outdir,paste("cna1ko_crz1ko_intersect",fileend,sep="_")))
    write.csv(setdiff(KOcna1Genes,KOcrz1Genes),file=file.path(outdir,paste("cna1ko_unique",fileend,sep="_")))
    write.csv(setdiff(KOcrz1Genes,KOcna1Genes),file=file.path(outdir,paste("crz1ko_unique",fileend,sep="_")))
    ## return(list("KOcna1Genes" = KOcna1Genes,"KOcrz1Genes"=KOcrz1Genes))
    return(ListOfGeneVecs)
}
##------- CRZ1 overexpression ------
Crz1OverexpressHeatmap = function(ddsHTSeq,outdir){
    cna1.id="CNAG_04796"
    crz1.id="CNAG_00156"
    crz1mat = counts(ddsHTSeq,normalized=TRUE)[c(cna1.id,crz1.id),]
    colnames(crz1mat) = with(colData(ddsHTSeq),paste(condition, temp, sep=" : "))
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    pdf(file.path(outdir,"crz1_overexpress.pdf"))
    heatmap.2(t(crz1mat), trace="none", col = rev(hmcol), margin=c(7, 7),cexCol=1)
    dev.off()
}

##------- GenesOfInterestHeatmap ------
GenesOfInterestHeatmap = function(genes.of.interest,ddsHTSeq,outfile){
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    # pdf(outfile)
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Raw counts
    png(paste(outfile,"_raw.png",sep=""),width = 8, height = 10, units = "in", res = 300)
    goimat = counts(ddsHTSeq,normalized=TRUE)[genes.of.interest,]
    colnames(goimat) = with(colData(ddsHTSeq),paste(condition, temp, sep=" : "))
    heatmap.2(goimat, trace="none", col = rev(hmcol), margin=c(7, 7),cexCol=1)
    dev.off()
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Regularized log transformation
    png(paste(outfile,"_rld.png",sep=""),width = 8, height = 10, units = "in", res = 300)
    rld <- rlogTransformation(ddsHTSeq, blind=FALSE)
    goimat = assay(rld)[genes.of.interest,]
    colnames(goimat) = with(colData(rld),paste(condition, temp, sep=" : "))
    heatmap.2(goimat, trace="none", col = rev(hmcol), margin=c(7, 7),cexCol=1)
    dev.off()
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Regularized log transformation
    png(paste(outfile,"_vsd.png",sep=""),width = 8, height = 10, units = "in", res = 300)
    # rld <- rlogTransformation(ddsHTSeq, blind=TRUE)
    # vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
    vsd <- varianceStabilizingTransformation(ddsHTSeq, blind=FALSE)
    goimat = assay(vsd)[genes.of.interest,]
    colnames(goimat) = with(colData(vsd),paste(condition, temp, sep=" : "))
    heatmap.2(goimat, trace="none", col = rev(hmcol), margin=c(7, 7),cexCol=1)
    dev.off()
}
##----- Heatmap of the sample-to-sample distances------
SampleSampleDistHeatmap = function(ddsHTSeq,outdir){
    ## select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30] 
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

    rld <- rlogTransformation(ddsHTSeq, blind=TRUE)
    distsRL <- dist(t(assay(rld)))
    ## A heatmap of this distance matrix gives us an overview over similarities and dissimilarities between samples (Figure 8):
    mat <- as.matrix(distsRL)
    rownames(mat) <- colnames(mat) <- with(colData(ddsHTSeq),
                                           paste(condition, temp, sep=" : "))
    pdf(sampleCompareHeatmap)
    heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(7, 7))
    dev.off()
    
    pdf(sampleComparePCA)
    sampleplt = plotPCA(rld, intgroup=c("condition", "temp"))
    print(sampleplt)
    dev.off()
}

##------- Normalized Counts Table ------
ExportNormedCounts = function(KOcna1Genes,KOcrz1Genes,ddsHTSeq,fileend,outdir){
    ## head(counts(ddsHTSeq,normalized=TRUE))
    ## with(colData(ddsHTSeq),paste(row.names(colData(ddsHTSeq)),paste(condition, temp, sep="_"),sep=":"))
    ddsHTSeq.counts = counts(ddsHTSeq,normalized=TRUE)
    coldat = colData(ddsHTSeq)
    colnames(ddsHTSeq.counts) = paste(row.names(coldat),paste(coldat$condition, coldat$temp, sep="_"),sep=":")

    # fileend=paste(fdrcutoff*100,"fdr_counts.csv",sep="")
    # fileend=paste(fdrcutoff*100,"fdr_", fccutoff,"fc_counts.csv",sep="")
    fileend=paste(fileend,"_counts.csv",sep="")


    write.csv(ddsHTSeq.counts,file=file.path(outdir,"cna1_crz1_allcounts.csv"))
    write.csv(ddsHTSeq.counts[intersect(KOcna1Genes,KOcrz1Genes),],
              file=file.path(outdir,paste("cna1ko_crz1ko_intersect",fileend,sep="_")))
    write.csv(ddsHTSeq.counts[setdiff(KOcna1Genes,KOcrz1Genes),],
              file=file.path(outdir,paste("cna1ko_unique",fileend,sep="_")))
    write.csv(ddsHTSeq.counts[setdiff(KOcrz1Genes,KOcna1Genes),],
              file=file.path(outdir,paste("crz1ko_unique",fileend,sep="_")))
}
##------- Results Table ------
## HERE >>>>>>
## ExportResults = function(KOcna1Genes,KOcrz1Genes,ddsHTSeq,fdrcutoff,fccutoff,outdir){
##     coldat = colData(ddsHTSeq)
##     colnames(ddsHTSeq.counts) = paste(row.names(coldat),paste(coldat$condition, coldat$temp, sep="_"),sep=":")

##     # fileend=paste(fdrcutoff*100,"fdr_counts.csv",sep="")
##     fileend=paste(fdrcutoff*100,"fdr_", fccutoff,"fc_counts.csv",sep="")

##     write.csv(ddsHTSeq.counts,file=file.path(outdir,"cna1_crz1_allcounts.csv"))
##     write.csv(ddsHTSeq.counts[intersect(KOcna1Genes,KOcrz1Genes),],
##               file=file.path(outdir,paste("cna1ko_crz1ko_intersect",fileend,sep="_")))
##     write.csv(ddsHTSeq.counts[setdiff(KOcna1Genes,KOcrz1Genes),],
##               file=file.path(outdir,paste("cna1ko_unique",fileend,sep="_")))
##     write.csv(ddsHTSeq.counts[setdiff(KOcrz1Genes,KOcna1Genes),],
##               file=file.path(outdir,paste("crz1ko_unique",fileend,sep="_")))
## }
## HERE <<<<<<
##------- Annotate Gene Lists ------
AnnotateGeneLists = function(KOcna1Genes,KOcrz1Genes,fileend,annotdir, outdir){
    ##------------ Load Annotation --------
    # annotdir = file.path(Sys.getenv("BREM"),"cneo_hybrid_rnaseq/info")
    
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

    # fileend=paste(fdrcutoff*100,"fdr_annot.csv",sep="")
    # fileend=paste(fdrcutoff*100,"fdr_", fccutoff,"fc_annot.csv",sep="")
    fileend=paste(fileend,"_annot.csv",sep="")

    write.csv(annot.df[intersect(KOcna1Genes,KOcrz1Genes),],
              file=file.path(outdir,paste("cna1ko_crz1ko_intersect",fileend,sep="_")))
    write.csv(annot.df[setdiff(KOcna1Genes,KOcrz1Genes),],
              file=file.path(outdir,paste("cna1ko_unique",fileend,sep="_")))
    write.csv(annot.df[setdiff(KOcrz1Genes,KOcna1Genes),],
              file=file.path(outdir,paste("crz1ko_unique",fileend,sep="_")))
}
##===========================================================================
##===========================================================================
## get.genes.counts = function(comparison,dds,fdr){
##     cur.res = results(dds,comparison)
##     # return(row.names(cur.res[which(cur.res$padj < fdr),]))
##     return(sum(cur.res$padj < fdr,na.rm = TRUE))
## }
## mapply(FUN=get.genes.counts,resultsNames(ddsHTSeq),MoreArgs = list(dds=ddsHTSeq,fdr=0.2))
## mapply(FUN=get.genes.counts,resultsNames(ddsHTSeq),MoreArgs = list(dds=ddsHTSeq,fdr=0.05))
##===========================================================================
##===========================================================================
Crz1OverexpressHeatmap(ddsHTSeq,outdir)

# fdrcutoff = 0.05
# fdrcutoff = 0.2
# fdrlist = c(0.2, 0.05)
# curfc = 0.01
curfc = 0 ##HERE
fdrlist = c(0.05)
countfilter = FALSE

for (curfdr in fdrlist){
    if (countfilter) {
        filtered="_filt"
    } else {
        filtered=""
    }
    ## fileend=paste(curfdr*100,"fdr_", curfc,"fc",filtered,".csv",sep="")
    fileend=paste(curfdr*100,"fdr_", curfc,"fc",filtered,sep="")
    
    fdg = FindDiffGenes(ddsHTSeq, outdir,fdrcutoff=curfdr,fccutoff=curfc)
    KOcna1Genes = fdg[["KO_cna1"]] ## fdg$KOcna1Genes
    KOcrz1Genes = fdg[["KO_crz1"]] ## fdg$KOcrz1Genes
    ExportNormedCounts(KOcna1Genes,KOcrz1Genes,ddsHTSeq,fileend,outdir)
    AnnotateGeneLists(KOcna1Genes,KOcrz1Genes,fileend,annotdir,outdir)
    ## ExportResults(KOcna1Genes,KOcrz1Genes,ddsHTSeq,curfdr,curfc,outdir) ##HERE
    cna1ko.crz1ko.intersect.genes = intersect(KOcna1Genes,KOcrz1Genes)
    cna1ko.unique.genes = setdiff(KOcna1Genes,KOcrz1Genes)
    crz1ko.unique = setdiff(KOcrz1Genes,KOcna1Genes)

    ## fileend=paste(curfdr*100,"fdr_", curfc,"fc.csv",sep="")
    GenesOfInterestHeatmap(c(cna1ko.unique.genes,cna1ko.crz1ko.intersect.genes,crz1ko.unique),
                           ddsHTSeq, 
                           outfile= file.path(outdir,paste("goi_heatmap_", fileend,sep="")))
    GenesOfInterestHeatmap(cna1ko.unique.genes,
                           ddsHTSeq, 
                           outfile= file.path(outdir,paste("cna1ko_unique_heatmap_", fileend,sep="")))
    GenesOfInterestHeatmap(cna1ko.crz1ko.intersect.genes,
                           ddsHTSeq, 
                           outfile= file.path(outdir,paste("cna1ko_crz1ko_intersect_heatmap_", fileend,sep="")))
    GenesOfInterestHeatmap(crz1ko.unique,
                           ddsHTSeq, 
                           outfile= file.path(outdir,paste("crz1ko_unique_heatmap_", fileend,sep="")))

}
## Do following last because it ??alters esitmates??
SampleSampleDistHeatmap(ddsHTSeq,outdir)

##------- Individual Comparisons ------
# WT_24C <-> KO_crz1_24C
# WT_24C <-> KO_cna1_24C
# KO_cna1_24C <->  KO_crz1_24C

# WT_37C <-> KO_crz1_37C
# WT_37C <-> KO_cna1_37C
# KO_cna1_37C <->  KO_crz1_37C

# KO_cna1_24C <->  KO_cna1_37C
# KO_crz1_24C <->  KO_crz1_37C
stop("need to do filtering for each comparison???? -> see FindDiffGenes")
stop("Output results tables in FindDiffGenes?")
