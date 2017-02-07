library("readxl")
library(dplyr)
library(magrittr)
library(stringr)

getwd()
expression_data_dir= "rewiring_analysis/fold_change_tables/"
list.files(expression_data_dir)

eve.excelfile = file.path(expression_data_dir,
                          "Gene datasets for A fumigatus and S cerevisiae Crz1 genes.xlsx")



#+ load_cn_data, echo=TRUE
cn.df = read_excel(eve.excelfile, sheet = 1, 
                     col_names = c("gene", "crz1D_logFC", "crz1D_padj", 
                                   "cna1D_logFC", "cna1D_padj"),
                     col_types = c("text", "numeric", "numeric", 
                                   "numeric", "numeric")) %>%
  filter(!is.na(crz1D_logFC))

##-----------------------------------------------------------------------------
#+ load_af_data, echo=TRUE
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
## "Afu2g05330" is a previous ID for "Afu2g05325", as per
## http://fungidb.org/fungidb/app/record/gene/Afu2g05325
head(af_crz_down.df)
tail(af_crz_down.df)

af_crz_up.df = read_excel(eve.excelfile, sheet = 5, col_names = c("gene"), skip=1)  %>%
  mutate(gene=str_trim(gene)) %>%
  filter(!gene %in% c("Afu6g12810","Afu6g11860","Afu1g00190","Afu2g07390"))
# Afu6g12810, Afu6g11860, Afu1g00190, and Afu2g07390 do not appear in FungiDB

head(af_crz_up.df)
tail(af_crz_up.df)

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
ortholog_data_dir= "rewiring_analysis/ortholog_tables/"
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

#+ eliminate_multiple_mappings
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


