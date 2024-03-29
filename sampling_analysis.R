#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Load Libraries, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
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

results_outfile = file.path("rewiring_analysis", "resampling_analysis.txt")
file.remove(results_outfile)
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
#+ Setup: Load C. neoformans - A. fumigatus ortholog data, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cn_af_ortho.df = read.delim(file.path(ortholog_data_dir,"h99_afumigatus_orthologs.tsv"),
                            stringsAsFactors = FALSE) %>%
  transmute(gene=X.Gene.ID.,
            ortholog=X.Input.Ortholog.s..,
            paralog_count=X.Paralog.count.,
            ortholog_count=X.Ortholog.count.)

cn_af_ortho.sep.df = cn_af_ortho.df %>% separate_rows(ortholog,sep=",")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Load C. neoformans - S. cerevisiae ortholog data, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cn_sc_ortho.df = read.delim(file.path(ortholog_data_dir,"h99_scerevisiae_orthologs.tsv"),
                            stringsAsFactors = FALSE) %>%
  transmute(gene=X.Gene.ID.,
            ortholog=X.Input.Ortholog.s..,
            paralog_count=X.Paralog.count.,
            ortholog_count=X.Ortholog.count.)

cn_sc_ortho.sep.df = cn_sc_ortho.df %>% separate_rows(ortholog,sep=",")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Load S. cerevisiae - A. fumigatus ortholog data, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sc_af_ortho.df = read.delim(file.path(ortholog_data_dir,"scerevisiae_afumigatus_orthologs.tsv"),
                            stringsAsFactors = FALSE) %>%
  transmute(gene=X.Gene.ID.,
            ortholog=X.Input.Ortholog.s..,
            paralog_count=X.Paralog.count.,
            ortholog_count=X.Ortholog.count.)

sc_af_ortho.sep.df = sc_af_ortho.df %>% separate_rows(ortholog,sep=",")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  
  overlap_genes = ortho_table %>% 
    filter(gene %in% a_genes) %>% 
    filter(ortholog %in% b_genes)
  
  num_overlap_genes = nrow(overlap_genes)
  
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
  
  return(overlap_genes)
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Run Analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Comparison 1.1, echo=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
organism_a = "A. fumigatus"
organism_b = "C. neoformans"
# #' `r kable(af_counts.df)`

cn_af_overlap.df = RunSamplingAnalysis(af_crz_genes, cn_genes, 
                                       cn_af_ortho.sep.df, 
                                       organism_a, organism_b,
                                       num_samples=1000,seed=1,
                                       outfile=results_outfile)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ### Calcineurin-dependent genes shared between `r organism_a` and `r organism_b`
#+ Comparison 1.2, echo=FALSE
cn_af_overlap.df %>% transmute(a=gene, b=ortholog) %>% kable


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Comparison 2.1, echo=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
organism_a = "A. fumigatus"
organism_b = "S. cerevisiae"
sc_af_overlap.df = RunSamplingAnalysis(af_crz_genes, sc_genes, 
                                       sc_af_ortho.sep.df, 
                                       organism_a, organism_b,
                                       num_samples=1000,seed=2,
                                       outfile=results_outfile)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ### Calcineurin-dependent genes shared between `r organism_a` and `r organism_b`
#+ Comparison 2.2, echo=FALSE
sc_af_overlap.df %>% transmute(a=gene, b=ortholog) %>% kable

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Comparison 3.1, echo=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
organism_a = "S. cerevisiae"
organism_b = "C. neoformans"
cn_sc_overlap.df = RunSamplingAnalysis(sc_genes, cn_genes, 
                                       cn_sc_ortho.sep.df,
                                       organism_a, organism_b,
                                       num_samples=1000,seed=3,
                                       outfile=results_outfile)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ### Calcineurin-dependent genes shared between `r organism_a` and `r organism_b`
#+ Comparison 3.2, echo=FALSE
cn_sc_overlap.df %>% transmute(a=gene, b=ortholog) %>% kable


#' ### Calcineurin-dependent genes shared among all three species
#+ 3-way comparison, echo=FALSE
inner_join(cn_af_overlap.df,cn_sc_overlap.df, by = "ortholog" ) %>% 
  transmute(a=gene.x, b=ortholog, c=gene.y) %>% kable

#' Note that there are two homologs of CNAG_01232	and YGL006W in A. fumigatus: Afu3g10690 and Afu7g01030
  

# inner_join(cn_af_overlap.df,sc_af_overlap.df, by = "gene" )

# inner_join(sc_af_overlap.df,cn_sc_overlap.df, by = c("ortholog","gene"))


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
