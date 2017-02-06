if (interactive()){
    basedir<<-file.path(Sys.getenv("CNA"),"rstudio")
} else {
    basedir<<-"."
}
##================================================================================
outdir=file.path(basedir,"results")
annotdir = file.path(basedir,"info")

suppressPackageStartupMessages(library("DESeq2",lib.loc="/Users/josh/Library/R/3.0/library"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("gplots"))
## writeLines(capture.output(sessionInfo()), file.path(outdir,"sessionInfo.txt"))
print(.libPaths())
print("==================================================")
print(sessionInfo())
