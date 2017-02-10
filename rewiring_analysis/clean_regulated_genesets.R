#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Load Libraries, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(readxl)
library(dplyr)
library(magrittr)
library(stringr)
library(knitr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Setup Paths, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getwd()
list.files("rewiring_analysis/")
expression_data_dir= "rewiring_analysis/fold_change_tables"
ortholog_data_dir= "rewiring_analysis/ortholog_tables"
eve.excelfile = file.path(expression_data_dir,
                          "Gene datasets for A fumigatus and S cerevisiae Crz1 genes.xlsx")
clean_reg_gene_dir= "rewiring_analysis/clean_reg_gene_lists"

dir.create(clean_reg_gene_dir)
clean_cn_reg_genes = file.path(clean_reg_gene_dir,"cneoformans_regulated_genes.csv")

clean_af_down_genes = file.path(clean_reg_gene_dir,"afumigatus_down_regulated_genes.csv")
clean_af_up_genes = file.path(clean_reg_gene_dir,"afumigatus_up_regulated_genes.csv")

clean_sc_ca_genes = file.path(clean_reg_gene_dir,"scerevisiae_ca_regulated_genes.csv")
clean_sc_ca_na_genes = file.path(clean_reg_gene_dir,"scerevisiae_ca_na_regulated_genes.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Load C. neoformans data from excel, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cn.df = read_excel(eve.excelfile, sheet = 1, skip=2,
                     col_names = c("gene", "crz1D_logFC", "crz1D_padj", 
                                   "cna1D_logFC", "cna1D_padj"),
                     col_types = c("text", "numeric", "numeric", 
                                   "numeric", "numeric")) %>%
  filter(!is.na(crz1D_logFC))

write.table(cn.df$gene,file=clean_cn_reg_genes,
            sep=",", row.names = FALSE,col.names = FALSE,quote=FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: A. fumigatus data from excel, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# af_ca_up.df = read_excel(eve.excelfile, sheet = 2, col_names = c("gene"), skip=1) %>%
#   mutate(gene=str_trim(gene))
# 
# af_ca_down.df = read_excel(eve.excelfile, sheet = 3, col_names = c("gene"), skip=1)  %>%
#   mutate(gene=str_trim(gene))

af_crz_down.df = read_excel(eve.excelfile, sheet = 4, col_names = c("gene"), skip=1) %>%
  mutate(gene=str_trim(gene)) %>%
  mutate(gene=str_replace(gene,"Afu2g05330","Afu2g05325"))

write.table(af_crz_down.df,file=clean_af_down_genes,
            sep=",", row.names = FALSE,col.names = FALSE,quote=FALSE)

af.unknown.genes = c("Afu6g12810","Afu6g11860","Afu1g00190","Afu2g07390")
af_crz_up.df = read_excel(eve.excelfile, sheet = 5, col_names = c("gene"), skip=1)  %>%
  mutate(gene=str_trim(gene)) %>%
  filter(!gene %in% af.unknown.genes)

write.table(af_crz_up.df,file=clean_af_up_genes,
            sep=",", row.names = FALSE,col.names = FALSE,quote=FALSE)


# af_counts.df = data_frame(
#   data_set= c(
#     "Up-regulated by Ca in WT",
#     "Down-regulated by Ca in WT",
#     "Up-regulated by Ca in crzA mutant",
#     "Down-regulated by Ca in crzA mutant"),
#   gene_counts = c(nrow(af_ca_up.df),
#                   nrow(af_ca_down.df),
#                   nrow(af_crz_up.df),
#                   nrow(af_crz_down.df)))

#' ## Loading *A. fumigatus* differentially expressed genes
#' Renamed 'Afu2g05330' to 'Afu2g05325', since  'Afu2g05330' is a previous ID 
#' for "Afu2g05325", according to the 
#' [FungiDB record for Afu2g05325](http://fungidb.org/fungidb/app/record/gene/Afu2g05325)
#' 
#' Removed the following genes from *A. fumigatus* differentially expressed gene set
#' because they do not appear in FungiDB, and all are listed as 
#' "hypothetical protein" in Soriani paper: `r af.unknown.genes`

#' Number of differentially expressed genes in each set 
#' (after removing 'unknown' genes):
#' `r kable(af_counts.df)`


##-----------------------------------------------------------------------------
#+ load_sc_data, echo=TRUE
sc_cna_ca.df = read_excel(eve.excelfile, sheet = 6, col_names = c("gene"), skip=2) %>%
  mutate(gene=str_extract(gene, pattern="\\w+")) %>% 
  filter(!is.na(gene)) %>%
  mutate(gene=str_replace(gene,"Y0","YO")) %>%
  mutate(gene=str_replace(gene,"YMB304C","YMR304C-A")) %>%
  filter(!gene %in% c("YNL043C","YMR007W","YPR197C","YMR304C-A","YDL011C","YDL172C","YPR170C")) # remove genes not in FungiDB

write.table(sc_cna_ca.df,file=clean_sc_ca_genes,
            sep=",", row.names = FALSE,col.names = FALSE,quote=FALSE)


sc_cna_ca_na.df = read_excel(eve.excelfile, sheet = 7, col_names = c("gene"), skip=2) %>%
  filter(!gene %in% c("YDR535C","YGL165C")) # remove genes not in FungiDB

write.table(sc_cna_ca_na.df,file=clean_sc_ca_na_genes,
            sep=",", row.names = FALSE,col.names = FALSE,quote=FALSE)

##-----------------------------------------------------------------------------
#+ # Validate Gene Lists

#+ ## Validate Af Gene Lists
af.allgenes = read.delim(file.path(ortholog_data_dir,"afumigatus_all_genes.tsv"),
                         stringsAsFactors=FALSE)
# all(af_ca_up.df$gene %in% af.allgenes[,1])
# all(af_ca_down.df$gene %in% af.allgenes[,1])
all(af_crz_down.df$gene %in% af.allgenes[,1])
all(af_crz_up.df$gene %in% af.allgenes[,1])
rm(af.allgenes)

#+ ## Validate Cn Gene Lists
cn.allgenes = read.delim(file.path(ortholog_data_dir,"h99_all_genes.tsv"),
                         stringsAsFactors=FALSE)
all(cn.df$gene %in% cn.allgenes[,1])
rm(cn.allgenes)

#+ ## Validate sc Gene Lists
sc.allgenes = read.delim(file.path(ortholog_data_dir,"scerevisiae_all_genes.tsv"),
                         stringsAsFactors=FALSE)
all(sc_cna_ca.df$gene %in% sc.allgenes[,1])
all(sc_cna_ca_na.df$gene %in% sc.allgenes[,1])
rm(sc.allgenes)

sessionInfo()
