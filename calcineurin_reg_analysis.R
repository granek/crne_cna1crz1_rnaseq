library("DESeq2")
library("RColorBrewer")
library("gplots")

##================================================================================
if (interactive()){
    basedir<<-file.path(Sys.getenv("CNA"),"rstudio")
} else {
    basedir<<-"."
}
##================================================================================
outdir=file.path(basedir,"results")
annotdir = file.path(basedir,"info")
countdir=file.path(basedir,"counts")
dir.create(outdir)
##================================================================================
suppressPackageStartupMessages(library("optparse"))

option_list <- list( 
    make_option("--table", default=file.path(basedir,"info","calcineurin_sample_table.csv"), 
                help = "Sample table file for count data [default: \"%default\"]"),
    make_option("--label", default="", 
                help = "Prefix LABEL to output file names"),
    make_option("--fc", type="numeric", default=2,
                help="Minimum fold change cutoff for genes"),
    make_option("--fdr", type="numeric", default=0.2,
                help="Maximum false discovery rate (FDR) cutoff for genes")
    )
opt <- parse_args(OptionParser(option_list=option_list))

# counttab.file=file.path(basedir,"info","calcineurin_sample_table_drop_bad_cnako.csv")
counttab.file=opt$table
outbase = file.path(outdir,opt$label)
##================================================================================
sampleData = read.csv(counttab.file,comment.char="#", colClasses=c("character","numeric","character","factor","factor"))
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



FindDiffGenes = function(ddsHTSeq,outbase, fdrcutoff=0.05, fccutoff=2,countfilter=FALSE){
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
        write.csv(as.data.frame(filt.res),file=paste(sep="",outbase,paste(sep="_",coeff, "results",fileend)))
        ListOfGeneVecs[[sample]] = row.names(filt.res)
    }

    KOcna1Genes = ListOfGeneVecs[["KO_cna1"]]
    KOcrz1Genes = ListOfGeneVecs[["KO_crz1"]]
    write.csv(KOcna1Genes,file=paste(sep="",outbase,paste("cna1ko_genes",fileend,sep="_")))
    write.csv(KOcrz1Genes,file=paste(sep="",outbase,paste("crz1ko_genes",fileend,sep="_")))
    write.csv(intersect(KOcna1Genes,KOcrz1Genes),file=paste(sep="",outbase,paste("cna1ko_crz1ko_intersect",fileend,sep="_")))
    write.csv(setdiff(KOcna1Genes,KOcrz1Genes),file=paste(sep="",outbase,paste("cna1ko_unique",fileend,sep="_")))
    write.csv(setdiff(KOcrz1Genes,KOcna1Genes),file=paste(sep="",outbase,paste("crz1ko_unique",fileend,sep="_")))
    ## return(list("KOcna1Genes" = KOcna1Genes,"KOcrz1Genes"=KOcrz1Genes))
    return(ListOfGeneVecs)
}
##------- CRZ1 overexpression ------
Crz1OverexpressHeatmap = function(ddsHTSeq,outbase,clustmethod="average"){
    cna1.id="CNAG_04796"
    crz1.id="CNAG_00156"
    crz1mat = counts(ddsHTSeq,normalized=TRUE)[c(cna1.id,crz1.id),]
    cols = colData(ddsHTSeq)
    ## colnames(crz1mat) = with(colData(ddsHTSeq),paste(condition, temp, sep=" : "))
    colnames(crz1mat) = paste(cols$condition, cols$temp, substr(rownames(cols),6,7), sep=" : ")
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    ## pdf(file.path(outdir,"crz1_overexpress.pdf"))
    png(paste(sep="",outbase,"crz1_overexpress.png"),width = 8, height = 10, units = "in", res = 300)
    heatmap.2(t(crz1mat), trace="none", col = rev(hmcol), margin=c(7, 7),cexCol=1,
              hclustfun=function(d,members=NULL){hclust(d,method=clustmethod,members)})
    dev.off()
}

##------- GenesOfInterestHeatmap ------
GenesOfInterestHeatmap = function(genes.of.interest,ddsHTSeq,outfile,clustmethod="average"){
    print(paste("NUMBER OF genes.of.interest:",length(genes.of.interest)))
    if (length(genes.of.interest) == 0){
        warning("Refusing to make a heatmap with ZERO genes.of.interest")
        return(FALSE)
    }
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    # pdf(outfile)
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Raw counts
    png(paste(outfile,"_raw.png",sep=""),width = 8, height = 10, units = "in", res = 300)
    goimat = counts(ddsHTSeq,normalized=TRUE)[genes.of.interest,]
    colnames(goimat) = with(colData(ddsHTSeq),paste(condition, temp, sep=" : "))
    heatmap.2(goimat, trace="none", col = rev(hmcol), margin=c(7, 7),cexRow=1,cexCol=1,
              hclustfun=function(d,members=NULL){hclust(d,method=clustmethod,members)})
    dev.off()
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Regularized log transformation
    png(paste(outfile,"_rld.png",sep=""),width = 8, height = 10, units = "in", res = 300)
    rld <- rlogTransformation(ddsHTSeq, blind=FALSE)
    goimat = assay(rld)[genes.of.interest,]
    colnames(goimat) = with(colData(rld),paste(condition, temp, sep=" : "))
    heatmap.2(goimat, trace="none", col = rev(hmcol), margin=c(7, 7),cexRow=1,cexCol=1,
              hclustfun=function(d,members=NULL){hclust(d,method=clustmethod,members)})
    dev.off()
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Regularized log transformation
    png(paste(outfile,"_vsd.png",sep=""),width = 8, height = 10, units = "in", res = 300)
    # rld <- rlogTransformation(ddsHTSeq, blind=TRUE)
    # vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
    vsd <- varianceStabilizingTransformation(ddsHTSeq, blind=FALSE)
    goimat = assay(vsd)[genes.of.interest,]
    colnames(goimat) = with(colData(vsd),paste(condition, temp, sep=" : "))
    heatmap.2(goimat, trace="none", col = rev(hmcol), margin=c(7, 7),cexRow=1,cexCol=1,
              hclustfun=function(d,members=NULL){hclust(d,method=clustmethod,members)})

    dev.off()
}
##----- Heatmap of the sample-to-sample distances------
SampleSampleDistHeatmap = function(ddsHTSeq,outbase,clustmethod="average"){
    ## select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30] 
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)


    ## samplecols = c("black","blue", "lightblue", "green","lightgreen") #brewer.pal(5, "Paired")
    samplecols = brewer.pal(5, "Paired")
    names(samplecols) = c("KO_crz1","KI_CRZ1","KO_cna1","KI_CNA1","WT")
    colorbar = samplecols[as.character(colData(ddsHTSeq)$condition)]
    print(colorbar)
    print(samplecols)
    print(colData(ddsHTSeq)$condition)
    rld <- rlogTransformation(ddsHTSeq, blind=TRUE)
    distsRL <- dist(t(assay(rld)))
    ## A heatmap of this distance matrix gives us an overview over similarities and dissimilarities between samples (Figure 8):
    mat <- as.matrix(distsRL)
    cols = colData(ddsHTSeq)
    ## colnames(crz1mat) = with(colData(ddsHTSeq),paste(condition, temp, sep=" : "))
    rownames(mat) <- colnames(mat) <- paste(cols$condition, cols$temp, substr(rownames(cols),6,7), sep=" : ")
    ## rownames(mat) <- colnames(mat) <- with(colData(ddsHTSeq),
    ##                                        paste(condition, temp, sep=" : "))
    ## sampleCompareHeatmap = file.path(outdir,paste("sample_compare_heatmap_",clustmethod,".pdf",sep=""))
    sampleCompareHeatmap = paste(sep="",outbase,paste("sample_compare_heatmap_",clustmethod,".png",sep=""))
    # pdf(sampleCompareHeatmap)
    png(sampleCompareHeatmap,width = 7, height = 7, units = "in", res = 300)
    heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(7, 7),
              ColSideColors=colorbar,
              hclustfun=function(d,members=NULL){hclust(d,method=clustmethod,members)})

    dev.off()
    ## sampleComparePCA = file.path(outdir,"sample_compare_pca.pdf")
    ## pdf(sampleComparePCA)
    sampleComparePCA = paste(sep="",outbase,"sample_compare_pca.png")
    png(sampleComparePCA,width = 8, height = 10, units = "in", res = 300)
    sampleplt = plotPCA(rld, intgroup=c("condition", "temp"))
    print(sampleplt)
    dev.off()
}

##------- Normalized Counts Table ------
## ExportNormedCounts = function(KOcna1Genes,KOcrz1Genes,ddsHTSeq,fileend,outdir){
ExportNormedCounts = function(GeneVec,outpath,ddsHTSeq){
    ## head(counts(ddsHTSeq,normalized=TRUE))
    ## with(colData(ddsHTSeq),paste(row.names(colData(ddsHTSeq)),paste(condition, temp, sep="_"),sep=":"))
    ddsHTSeq.counts = counts(ddsHTSeq,normalized=TRUE)
    coldat = colData(ddsHTSeq)
    colnames(ddsHTSeq.counts) = paste(row.names(coldat),paste(coldat$condition, coldat$temp, sep="_"),sep=":")

    write.csv(ddsHTSeq.counts[GeneVec,],file=outpath)
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
AnnotateGeneLists = function(GeneVec,outpath,annotdir){
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

    write.csv(annot.df[GeneVec,],file=outpath)
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
Crz1OverexpressHeatmap(ddsHTSeq,outbase)

# fdrcutoff = 0.05
# fdrcutoff = 0.2
# fdrlist = c(0.2, 0.05)
# curfc = 0.01
# curfc = 0 ##HERE
## fclist = c(2,0)
## fdrlist = c(0.2, 0.05)
## filtlist = c(TRUE,FALSE)
## ##------------------------------------------------------------
## fclist = c(2)
## fdrlist = c(0.2)
## ##------------------------------------------------------------
##------------------------------------------------------------
fclist = c(opt$fc)
fdrlist = c(opt$fdr)
##------------------------------------------------------------
filtlist = c(TRUE)
for (curfc in fclist) {
    for (countfilter in filtlist){
        for (curfdr in fdrlist){
            if (countfilter) {
                filtered="_filt"
            } else {
                filtered=""
            }
            fileend=paste(curfdr*100,"fdr_", curfc,"fc",filtered,sep="")
            
            fdg = FindDiffGenes(ddsHTSeq, outbase,fdrcutoff=curfdr,fccutoff=curfc,countfilter=countfilter)
            KOcna1Genes = fdg[["KO_cna1"]] ## fdg$KOcna1Genes
            KOcrz1Genes = fdg[["KO_crz1"]] ## fdg$KOcrz1Genes

            # ExportNormedCounts(KOcna1Genes,KOcrz1Genes,ddsHTSeq,fileend,outdir)
            ## AnnotateGeneLists(KOcna1Genes,KOcrz1Genes,fileend,annotdir,outdir)
            ## ExportResults(KOcna1Genes,KOcrz1Genes,ddsHTSeq,curfdr,curfc,outdir) ##HERE
            ##-----------------------------------------------------------------
            ListOfGeneSets = list()
            ListOfGeneSets[["cna1ko_crz1ko_intersect"]] = intersect(KOcna1Genes,KOcrz1Genes)
            ListOfGeneSets[["cna1ko_unique"]] = setdiff(KOcna1Genes,KOcrz1Genes)
            ListOfGeneSets[["crz1ko_unique"]] = setdiff(KOcrz1Genes,KOcna1Genes)
            GenesOfInterestHeatmap(stack(ListOfGeneSets)$values,ddsHTSeq, 
                                   outfile= paste(sep="",outbase,paste("goi_heatmap_", fileend,sep="")))
            ListOfGeneSets[["KI_CRZ1"]] = fdg[["KI_CRZ1"]]

            ExportNormedCounts(rownames(ddsHTSeq),
                               paste(sep="",outbase,"cna1_crz1_allcounts.csv"),
                               ddsHTSeq)
            for (name in names(ListOfGeneSets)) {
                GenesOfInterestHeatmap(ListOfGeneSets[[name]],ddsHTSeq, 
                                        outfile= paste(sep="",outbase,paste(name,"_heatmap_", fileend,sep="")))
                AnnotateGeneLists(ListOfGeneSets[[name]],
                                  paste(sep="",outbase,paste(name,fileend,"annot.csv",sep="_")),
                                  annotdir)
                ## ExportNormedCounts(KOcna1Genes,KOcrz1Genes,ddsHTSeq,fileend,outdir)
                ## fileend=paste(fileend,"_counts.csv",sep="")
                ExportNormedCounts(ListOfGeneSets[[name]],
                                  paste(sep="",outbase,paste(name,fileend,"counts.csv",sep="_")),
                                  ddsHTSeq)
            }
        }
    }   
}
## Do following last because it ??alters esitmates??
for (clust in c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")) {
    SampleSampleDistHeatmap(ddsHTSeq,outbase,clustmethod=clust)
}
##------- Individual Comparisons ------
# WT_24C <-> KO_crz1_24C
# WT_24C <-> KO_cna1_24C
# KO_cna1_24C <->  KO_crz1_24C

# WT_37C <-> KO_crz1_37C
# WT_37C <-> KO_cna1_37C
# KO_cna1_37C <->  KO_crz1_37C

# KO_cna1_24C <->  KO_cna1_37C
# KO_crz1_24C <->  KO_crz1_37C
## stop("need to do filtering for each comparison???? -> see FindDiffGenes")
warning("Output results tables in FindDiffGenes?")
## stop("Use different clustering method")


test = function(x){
    if (x) {
        warning("blah1")
        return("now")
    }
    print("NOT x")
    warning("blah2")
}
