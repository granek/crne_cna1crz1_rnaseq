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

results_outfile = file.path("rewiring_analysis", "resampling_analysis.txt")
file.remove(results_outfile)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Load C. neoformans data from excel, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cn.df = read_excel(eve.excelfile, sheet = 1, 
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

##-----------------------------------------------------------------------------
list.files("rewiring_analysis/ortholog_tables/")

cn_af_ortho.df = read.delim(file.path(ortholog_data_dir,"h99_afumigatus_orthologs.tsv"),
                            stringsAsFactors = FALSE) %>%
                              transmute(gene=X.Gene.ID.,
                                        ortholog=X.Input.Ortholog.s..,
                                        paralog_count=X.Paralog.count.,
                                        ortholog_count=X.Ortholog.count.)

cn_af.nodups = cn_af_ortho.df %>%
  filter(paralog_count == 0) %>%
  filter(!str_detect(ortholog, ","))
##-----------------------------------------------------------------------------
#+ explore Af-Cn overlapping genes
af_ca.genes = c(af_ca_up.df$gene, af_ca_down.df$gene)
# cn_af.nodups %>% 
#   filter(gene %in% af_ca.genes) %>%
#   filter(ortholog %in% cn.df$gene)
# 
# cn_af.nodups %>% 
#   filter(ortholog %in% cn.df$gene)
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
sc_genes = c(sc_cna_ca.df$gene, sc_cna_ca_na.df$gene)

# cn_sc.nodups %>% 
#   filter(gene %in% sc_genes) %>%
#   filter(ortholog %in% cn.df$gene)
# 
# cn_af.nodups %>% 
#   filter(ortholog %in% cn.df$gene)

##-----------------------------------------------------------------------------
#+ Separate Cn-Sc paralogs into separate rows
cn_sc_ortho.sep.df = cn_sc_ortho.df %>% separate_rows(ortholog,sep=",")

##-----------------------------------------------------------------------------
#+ explore Cn-Sc overlapping genes
# cn_sc_ortho.sep.df %>% 
#   filter(gene %in% sc_genes) %>%
#   filter(ortholog %in% cn.df$gene)
# 
# cn_sc_ortho.sep.df %>% 
#   filter(gene %in% sc_genes) %>% dim
# 
# cn_sc_ortho.sep.df %>% 
#   filter(ortholog %in% cn.df$gene) %>% dim

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
#+ Separate Af-Cn paralogs into separate rows
cn_af_ortho.sep.df = cn_af_ortho.df %>% separate_rows(ortholog,sep=",")
# dim(cn_af_ortho.df)
# dim(cn_af_ortho.sep.df)
# 
# cn_af_ortho.df %>% filter(str_detect(ortholog, ",")) %>% head
# 
# cn_af_ortho.sep.df %>% filter(ortholog == "CNAG_00688")
##-----------------------------------------------------------------------------
#+ explore Af-Cn overlapping genes redux
af_crz.genes = c(af_crz_up.df$gene, af_crz_down.df$gene)
# cn_af_ortho.sep.df %>% 
#   filter(gene %in% af_crz.genes) %>%
#   filter(ortholog %in% cn.df$gene)
# 
# cn_af.nodups %>% 
#   filter(gene %in% af_crz.genes) %>%
#   filter(ortholog %in% cn.df$gene)


RunSamplingAnalysis = function(a_genes, b_genes, ortho_table, 
                               a_species, b_species,num_samples=100,
                               seed=NULL,outfile=""){
  if(!is.null(seed)){
    set.seed(seed)
  }
  a_num_genes = length(a_genes)
  a_num_genes_with_ortholog = ortho_table %>% 
    filter(gene %in% a_genes) %>% 
    nrow
  
  b_num_genes = length(b_genes)
  b_num_genes_with_ortholog = ortho_table %>% 
    filter(ortholog %in% b_genes) %>% 
    nrow
  
  num_overlap_genes = ortho_table %>% 
    filter(gene %in% a_genes) %>% 
    filter(ortholog %in% b_genes) %>% 
    nrow
  
  a_orthologs = unique(ortho_table$gene)
  b_orthologs = unique(ortho_table$ortholog)
  map_table = ortho_table %>%
    transmute(a = gene,
              b = ortholog)
  
  
  cat(paste("Orthologs between", a_species, "and", b_species,
            ":", nrow(ortho_table), "\n"), fill=TRUE,file=outfile,append=TRUE)
  cat(paste("Differentially regulated genes in", a_species, ":", 
            a_num_genes), fill=TRUE,file=outfile,append=TRUE)
  cat(paste("Number of regulated genes in", a_species, "with ortholog in",
            b_species, ":",
            a_num_genes_with_ortholog, "\n"), fill=TRUE,file=outfile,append=TRUE)

  cat(paste("Differentially regulated genes in", b_species, ":", 
            b_num_genes), fill=TRUE,file=outfile,append=TRUE)
  cat(paste("Number of regulated genes in", b_species, "with ortholog in",
            a_species, ":",
            b_num_genes_with_ortholog,"\n"), fill=TRUE,file=outfile,append=TRUE)
  cat(paste("Number of shared regulated genes between", a_species, "and",
            b_species, ":",
            num_overlap_genes), fill=TRUE,file=outfile,append=TRUE)
  
  y = RepWrap(num_samples, a_genes = a_orthologs, b_genes = b_orthologs,
              a_gene_n = a_num_genes_with_ortholog,
              b_gene_n = b_num_genes_with_ortholog,
              map_table = map_table)
  sample_prob = mean(y >= num_overlap_genes)

  cat(paste("Probability of finding number of shared regulated genes in", a_species, "and",
            b_species, ":",
            sample_prob), fill=TRUE,file=outfile,append=TRUE)
  cat(paste0(rep("-", 60),collapse=""),fill=TRUE,file=outfile,append=TRUE)
}

SampleFunc = function(a_genes, b_genes, a_gene_n, b_gene_n, map_table){
  a_samp = sample(a_genes,a_gene_n)
  b_samp = sample(b_genes,b_gene_n)
  samp_overlap = map_table %>%
    filter(a %in% a_samp) %>%
    filter(b %in% b_samp) %>%
    nrow
  return(samp_overlap)
}
RepWrap = function(n, a_genes, b_genes, a_gene_n, b_gene_n, map_table) {
  replicate(n, SampleFunc(a_genes=a_genes, b_genes=b_genes, 
                          a_gene_n=a_gene_n, b_gene_n=b_gene_n, 
                          map_table=map_table))
}

RunSamplingAnalysis(af_crz.genes, cn.df$gene, cn_af_ortho.sep.df, 
                    "A. fumigatus", "C. neoformans",num_samples=1000,seed=1,
                    outfile=results_outfile)

# RunSamplingAnalysis(af_crz.genes, cn.df$gene, cn_af.nodups, 
#                     "A. fumigatus", "C. neoformans",num_samples=1000,seed=2)

RunSamplingAnalysis(sc_genes, cn.df$gene, cn_sc_ortho.sep.df,
                    "S. cerevisiae", "C. neoformans",num_samples=1000,seed=3,
                    outfile=results_outfile)

# RunSamplingAnalysis(sc_genes, cn.df$gene, cn_sc.nodups,
#                     "S. cerevisiae", "C. neoformans",num_samples=1000,seed=4)
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
