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
##================================================================================
library("DESeq2")
## library("RColorBrewer")
## library("gplots")
writeLines(capture.output(sessionInfo()), file.path(outdir,"sessionInfo.txt"))
##================================================================================
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
ddsHTSeq <- DESeq(ddsHTSeq,betaPrior = FALSE)
res <- results(ddsHTSeq)

head(res)
print("----------------------------------------")
print("----------------------------------------\nresultsNames")
print("----------------------------------------")
print(resultsNames(ddsHTSeq))
resultsNames(ddsHTSeq)[3:6]



FindDiffGenes = function(ddsHTSeq,outbase, fdrcutoff=0.05, fccutoff=2,countfilter=FALSE){
    ## stop("need to do filtering for each comparison????")
    log2fc = log2(fccutoff)
    # for(var in seq) expr
    ListOfGeneVecs = list()
    for(sample in c("KI_CNA1","KI_CRZ1","KO_cna1","KO_crz1")) {
        coeff = paste("condition",sample, "vs_WT",sep="_")
        outfile = paste(sep="_",coeff, "results",fileend)
        ## if (packageVersion("DESeq2")=="1.2.10"){
        ##     coeff = paste("condition",sample, "vs_WT",sep="_")
        ##     outfile = paste(sep="_",coeff, "results",fileend)
        ## }else{
        ##     coeff = paste("condition",sample,sep="")
        ##     outfile = paste("condition",sample, "vs_WT","results",fileend,sep="_")
        ## }
        print("coeff")
        print(coeff)
        cur.res = results(ddsHTSeq,name=coeff)
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
        outfile = paste("condition",sample, "vs_WT","results",fileend,sep="_")
        write.csv(as.data.frame(filt.res),file=paste(sep="",outbase,outfile))
        ListOfGeneVecs[[sample]] = row.names(filt.res)
    }

    KOcna1Genes = ListOfGeneVecs[["KO_cna1"]]
    KOcrz1Genes = ListOfGeneVecs[["KO_crz1"]]
    write.csv(KOcna1Genes,file=paste(sep="",outbase,paste("cna1ko_genes",fileend,sep="_")))
    write.csv(KOcrz1Genes,file=paste(sep="",outbase,paste("crz1ko_genes",fileend,sep="_")))
    ## write.csv(intersect(KOcna1Genes,KOcrz1Genes),file=paste(sep="",outbase,paste("cna1ko_crz1ko_intersect",fileend,sep="_")))
    ## write.csv(setdiff(KOcna1Genes,KOcrz1Genes),file=paste(sep="",outbase,paste("cna1ko_unique",fileend,sep="_")))
    ## write.csv(setdiff(KOcrz1Genes,KOcna1Genes),file=paste(sep="",outbase,paste("crz1ko_unique",fileend,sep="_")))
    ## return(list("KOcna1Genes" = KOcna1Genes,"KOcrz1Genes"=KOcrz1Genes))
    return(ListOfGeneVecs)
}

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
        }
    }   
}
