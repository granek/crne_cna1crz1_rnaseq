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
ortholog_data_dir= "rewiring_analysis/ortholog_tables"
clean_reg_gene_dir= "rewiring_analysis/clean_reg_gene_lists"

clean_cn_reg_genes = file.path(clean_reg_gene_dir,"cneoformans_regulated_genes.csv")

clean_af_down_genes = file.path(clean_reg_gene_dir,"afumigatus_down_regulated_genes.csv")
clean_af_up_genes = file.path(clean_reg_gene_dir,"afumigatus_up_regulated_genes.csv")

clean_sc_ca_genes = file.path(clean_reg_gene_dir,"scerevisiae_ca_regulated_genes.csv")
clean_sc_ca_na_genes = file.path(clean_reg_gene_dir,"scerevisiae_ca_na_regulated_genes.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Load Regulated Gene Lists
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
af_ca_up = scan(clean_af_up_genes, what=character())
af_ca_down = scan(clean_af_down_genes, what=character())
af_crz_genes = c(af_ca_up, af_ca_down)

sc_cna_ca = scan(clean_sc_ca_genes, what=character())
sc_cna_ca_na = scan(clean_sc_ca_na_genes, what=character())
sc_genes = c(sc_cna_ca, sc_cna_ca_na)

cn_genes = scan(clean_cn_reg_genes, what=character())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Load C. neoformans data from excel, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##-----------------------------------------------------------------------------
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
# cn_af.nodups %>% 
#   filter(gene %in% af_ca.genes) %>%
#   filter(ortholog %in% cn_genes)
# 
# cn_af.nodups %>% 
#   filter(ortholog %in% cn_genes)
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

# cn_sc.nodups %>% 
#   filter(gene %in% sc_genes) %>%
#   filter(ortholog %in% cn_genes)
# 
# cn_af.nodups %>% 
#   filter(ortholog %in% cn_genes)

##-----------------------------------------------------------------------------
#+ Separate Cn-Sc paralogs into separate rows
cn_sc_ortho.sep.df = cn_sc_ortho.df %>% separate_rows(ortholog,sep=",")

##-----------------------------------------------------------------------------
#+ explore Cn-Sc overlapping genes
# cn_sc_ortho.sep.df %>% 
#   filter(gene %in% sc_genes) %>%
#   filter(ortholog %in% cn_genes)
# 
# cn_sc_ortho.sep.df %>% 
#   filter(gene %in% sc_genes) %>% dim
# 
# cn_sc_ortho.sep.df %>% 
#   filter(ortholog %in% cn_genes) %>% dim

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
# cn_af_ortho.sep.df %>% 
#   filter(gene %in% af_crz_genes) %>%
#   filter(ortholog %in% cn_genes)
# 
# cn_af.nodups %>% 
#   filter(gene %in% af_crz_genes) %>%
#   filter(ortholog %in% cn_genes)


RunSamplingAnalysis = function(a_genes, b_genes, ortho_table, 
                               a_species, b_species,num_samples=100,seed=NULL){
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
            ":", nrow(ortho_table), "\n"), fill=TRUE)
  cat(paste("Differentially regulated genes in", a_species, ":", 
            a_num_genes), fill=TRUE)
  cat(paste("Number of regulated genes in", a_species, "with ortholog in",
            b_species, ":",
            a_num_genes_with_ortholog, "\n"), fill=TRUE)

  cat(paste("Differentially regulated genes in", b_species, ":", 
            b_num_genes), fill=TRUE)
  cat(paste("Number of regulated genes in", b_species, "with ortholog in",
            a_species, ":",
            b_num_genes_with_ortholog,"\n"), fill=TRUE)
  cat(paste("Number of shared regulated genes between", a_species, "and",
            b_species, ":",
            num_overlap_genes), fill=TRUE)
  
  y = RepWrap(num_samples, a_genes = a_orthologs, b_genes = b_orthologs,
              a_gene_n = a_num_genes_with_ortholog,
              b_gene_n = b_num_genes_with_ortholog,
              map_table = map_table)
  sample_prob = mean(y >= num_overlap_genes)

  cat(paste("Probability of finding number of shared regulated genes in", a_species, "and",
            b_species, ":",
            sample_prob), fill=TRUE)
  cat(paste0(rep("-", 60),collapse=""),fill=TRUE)
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

RunSamplingAnalysis(af_crz_genes, cn_genes, cn_af_ortho.sep.df, 
                    "A. fumigatus", "C. neoformans",num_samples=1000,seed=1)

# RunSamplingAnalysis(af_crz_genes, cn_genes, cn_af.nodups, 
#                     "A. fumigatus", "C. neoformans",num_samples=1000,seed=2)

RunSamplingAnalysis(sc_genes, cn_genes, cn_sc_ortho.sep.df,
                    "S. cerevisiae", "C. neoformans",num_samples=1000,seed=3)

# RunSamplingAnalysis(sc_genes, cn_genes, cn_sc.nodups,
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
