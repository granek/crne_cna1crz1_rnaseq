#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Load Libraries, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(readxl)
library(dplyr)
library(magrittr)
library(stringr)
library(tidyr)
library(knitr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Setup Paths, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getwd()
expression_data_dir= "fold_change_tables/"
eve.excelfile = file.path(expression_data_dir,
                          "Gene datasets for A fumigatus and S cerevisiae Crz1 genes.xlsx")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Load C. neoformans data from excel, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cn.df = read_excel(eve.excelfile, sheet = 1, 
                     col_names = c("gene", "crz1D_logFC", "crz1D_padj", 
                                   "cna1D_logFC", "cna1D_padj"),
                     col_types = c("text", "numeric", "numeric", 
                                   "numeric", "numeric")) %>%
  filter(!is.na(crz1D_logFC))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: A. fumigatus data from excel, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
af_ca_up.df = read_excel(eve.excelfile, sheet = 2, col_names = c("gene"), skip=1) %>%
  mutate(gene=str_trim(gene))
head(af_ca_up.df)
tail(af_ca_up.df)

af_ca_down.df = read_excel(eve.excelfile, sheet = 3, col_names = c("gene"), skip=1)  %>%
  mutate(gene=str_trim(gene))
head(af_ca_down.df)
tail(af_ca_down.df)

af_crz_down.df = read_excel(eve.excelfile, sheet = 4, col_names = c("gene"), skip=1) %>%
  mutate(gene=str_trim(gene)) %>%
  mutate(gene=str_replace(gene,"Afu2g05330","Afu2g05325")) 

head(af_crz_down.df)
tail(af_crz_down.df)

af.unknown.genes = c("Afu6g12810","Afu6g11860","Afu1g00190","Afu2g07390")
af_crz_up.df = read_excel(eve.excelfile, sheet = 5, col_names = c("gene"), skip=1)  %>%
  mutate(gene=str_trim(gene)) %>%
  filter(!gene %in% af.unknown.genes)

af_counts.df = data_frame(
  data_set= c(
    "Up-regulated by Ca in WT",
    "Down-regulated by Ca in WT",
    "Up-regulated by Ca in crzA mutant",
    "Down-regulated by Ca in crzA mutant"),
  gene_counts = c(nrow(af_ca_up.df),
                  nrow(af_ca_down.df),
                  nrow(af_crz_up.df),
                  nrow(af_crz_down.df)))

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

head(sc_cna_ca.df)
tail(sc_cna_ca.df)

sc_cna_ca_na.df = read_excel(eve.excelfile, sheet = 7, col_names = c("gene"), skip=2) %>%
  filter(!gene %in% c("YDR535C","YGL165C")) # remove genes not in FungiDB

head(sc_cna_ca_na.df)
tail(sc_cna_ca_na.df)

##-----------------------------------------------------------------------------
#+ # Validate Gene Lists
ortholog_data_dir= "ortholog_tables/"
list.files(ortholog_data_dir)

#+ ## Validate Af Gene Lists
af.allgenes = read.delim(file.path(ortholog_data_dir,"afumigatus_all_genes.tsv"),
                         stringsAsFactors=FALSE)
all(af_ca_up.df$gene %in% af.allgenes[,1])
all(af_ca_down.df$gene %in% af.allgenes[,1])
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

# head(sc_cna_ca.df)
# head(sc_cna_ca_na.df)
# all(sc_cna_ca.df$gene %in% sc.allgenes[,1])
# sc_cna_ca.df$gene %in% sc.allgenes[,1]
# sc_cna_ca.df[which(!sc_cna_ca.df$gene %in% sc.allgenes[,1]),]
# which(!sc_cna_ca.df$gene %in% sc.allgenes[,1])


rm(sc.allgenes)

# af_crz_up.df$gene[which(not(af_crz_up.df$gene %in% af.allgenes[,1]))]
# 
# af_ca_down.df
# 
# head(af.allgenes[,1])
# head(af_ca_up.df$gene)
# head(af_ca_up.df)
# af_crz_up.df[which(af_crz_up.df$gene %in% af.allgenes[,1]),]
# dim(af_crz_up.df)

##-----------------------------------------------------------------------------
list.files("rewiring_analysis/ortholog_tables/")

cn_af_ortho.df = read.delim(file.path(ortholog_data_dir,"h99_afumigatus_orthologs.tsv"),
                            stringsAsFactors = FALSE) %>%
                              transmute(gene=X.Gene.ID.,
                                        ortholog=X.Input.Ortholog.s..,
                                        paralog_count=X.Paralog.count.,
                                        ortholog_count=X.Ortholog.count.)

#+ Af-Cn eliminate_multiple_mappings
# head(cn_af_ortho.df)
# sum(duplicated(cn_af_ortho.df$ortholog))
# head(cn_af_ortho.df[duplicated(cn_af_ortho.df$ortholog),])

# head(cn_af_ortho.df$ortholog)
# all_cn_orthos = unlist(str_split(cn_af_ortho.df$ortholog,","))
# all_cn_orthos.dups = all_cn_orthos[duplicated(all_cn_orthos)]
# ortho_pattern = paste(all_cn_orthos.dups,collapse = "|")
# 
# # cn_af.nodups = cn_af_ortho.df %>% 
# #   mutate(dups = str_detect(ortholog, ortho_pattern)) %>%
# #   mutate(comma = str_detect(ortholog, ",")) %>%
# #   filter(dups) %>%
# #   filter(comma)
# # range(cn_af.nodups$paralog_count)
# 
# cn_af.nodups1 = cn_af_ortho.df %>% 
#   filter(!str_detect(ortholog, ortho_pattern)) %>%
#   filter(!str_detect(ortholog, ","))
#   
# range(cn_af.nodups$paralog_count)
# 
# cn_af.nodups2 = cn_af_ortho.df %>% 
#   filter(paralog_count == 0) %>%
#   filter(!str_detect(ortholog, ","))
# 
# cn_af.nodups3 = cn_af_ortho.df %>% 
#   filter(!str_detect(ortholog, ortho_pattern))
# 
# 
# all.equal(cn_af.nodups3, cn_af.nodups2)
# 
# length(all_cn_orthos[duplicated(all_cn_orthos)])

cn_af.nodups = cn_af_ortho.df %>%
  filter(paralog_count == 0) %>%
  filter(!str_detect(ortholog, ","))
##-----------------------------------------------------------------------------
#+ explore Af-Cn overlapping genes
af_ca.genes = c(af_ca_up.df$gene, af_ca_down.df$gene)
cn_af.nodups %>% 
  filter(gene %in% af_ca.genes) %>%
  filter(ortholog %in% cn.df$gene)

cn_af.nodups %>% 
  filter(ortholog %in% cn.df$gene)


head(cn.df )
##-----------------------------------------------------------------------------
#+ Sc-Cn eliminate_multiple_mappings
list.files("rewiring_analysis/ortholog_tables/")

cn_sc_ortho.df = read.delim(file.path(ortholog_data_dir,"h99_scerevisiae_orthologs.tsv"),
                            stringsAsFactors = FALSE) %>%
  transmute(gene=X.Gene.ID.,
            ortholog=X.Input.Ortholog.s..,
            paralog_count=X.Paralog.count.,
            ortholog_count=X.Ortholog.count.)

cn_sc.nodups = cn_sc_ortho.df %>%
  filter(paralog_count == 0) %>%
  filter(!str_detect(ortholog, ","))

##-----------------------------------------------------------------------------
#+ explore Sc-Cn overlapping genes
sc_genes = c(sc_cna_ca.df$gene, sc_cna_ca_na.df)
# sc_cna_ca_na.df$gene
# dim(sc_cna_ca_na.df)
# dim(sc_cna_ca.df)

cn_sc.nodups %>% 
  filter(gene %in% sc_genes) %>%
  filter(ortholog %in% cn.df$gene)

cn_af.nodups %>% 
  filter(ortholog %in% cn.df$gene)


head(cn.df )

##-----------------------------------------------------------------------------
#+ Separate Cn-Sc paralogs into separate rows
head(cn_sc_ortho.df)
cn_sc_ortho.sep.df = cn_sc_ortho.df %>% separate_rows(ortholog,sep=",")
dim(cn_sc_ortho.df)
dim(cn_sc_ortho.sep.df)

cn_sc_ortho.df %>% filter(str_detect(ortholog, ",")) %>% head

cn_sc_ortho.sep.df %>% filter(gene == "YAL005C")

##-----------------------------------------------------------------------------
#+ explore Cn-Sc overlapping genes
cn_sc_ortho.sep.df %>% 
  filter(gene %in% sc_genes) %>%
  filter(ortholog %in% cn.df$gene)

cn_sc_ortho.sep.df %>% 
  filter(gene %in% sc_genes) %>% dim

cn_sc_ortho.sep.df %>% 
  filter(ortholog %in% cn.df$gene) %>% dim

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
#+ Separate Af-Cn paralogs into separate rows
cn_af_ortho.sep.df = cn_af_ortho.df %>% separate_rows(ortholog,sep=",")
dim(cn_af_ortho.df)
dim(cn_af_ortho.sep.df)

cn_af_ortho.df %>% filter(str_detect(ortholog, ",")) %>% head

cn_af_ortho.sep.df %>% filter(ortholog == "CNAG_00688")
##-----------------------------------------------------------------------------
#+ explore Af-Cn overlapping genes redux
af_crz.genes = c(af_crz_up.df$gene, af_crz_down.df$gene)
cn_af_ortho.sep.df %>% 
  filter(gene %in% af_crz.genes) %>%
  filter(ortholog %in% cn.df$gene)

cn_af.nodups %>% 
  filter(gene %in% af_crz.genes) %>%
  filter(ortholog %in% cn.df$gene)


RunSamplingAnalysis = function(a_genes, b_genes, ortho_table, a_species, b_species){
  num_a_genes = length(a_genes)
  num_a_genes_with_ortholog = ortho_table %>% 
    filter(gene %in% a_genes) %>% 
    nrow
  
  num_b_genes = length(b_genes)
  num_b_genes_with_ortholog = ortho_table %>% 
    filter(ortholog %in% b_genes) %>% 
    nrow
  
  num_overlap_genes = ortho_table %>% 
    filter(gene %in% a_genes) %>% 
    filter(ortholog %in% b_genes) %>% 
    nrow
  
  
  cat(paste("Orthologs between", a_species, "and", b_species,
      ":", nrow(ortho_table)), fill=TRUE)
  cat(paste("Differentially regulated genes in", a_species, ":", 
            num_a_genes), fill=TRUE)
  cat(paste("Number of regulated genes in", a_species, "with ortholog in",
            b_species, ":",
            num_a_genes_with_ortholog), fill=TRUE)

  cat(paste("Differentially regulated genes in", b_species, ":", 
            num_b_genes), fill=TRUE)
  cat(paste("Number of regulated genes in", b_species, "with ortholog in",
            a_species, ":",
            num_b_genes_with_ortholog), fill=TRUE)
  cat(paste("Number of common regulated genes between", a_species, "and",
            a_species, ":",
            num_overlap_genes), fill=TRUE)
}

RunSamplingAnalysis(af_crz.genes, cn.df$gene, cn_af_ortho.sep.df, 
                    "A. fumigatus", "C. neoformans")
#'******************************************************************************
#' # Further Analyses
#+ Todo List, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 1. One
#'     + sub1
#'     + sub2
#' 1. Two

#--------------------------------------------------
#'****************************************************************************
#' # Session Info
#--------------------------------------------------
#+ Session Info, echo=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sessionInfo()
