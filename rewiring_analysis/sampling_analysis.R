library("readxl")
library(dplyr)
library(magrittr)

getwd()
data.dir= "rewiring_analysis/fold_change_tables/"
list.files(data.dir)

eve.excelfile = file.path(data.dir,
                          "Gene datasets for A fumigatus and S cerevisiae Crz1 genes.xlsx")



#+ load_sheet1_crne, echo=TRUE
cn.df = read_excel(eve.excelfile, sheet = 1, 
                     col_names = c("gene", "crz1D_logFC", "crz1D_padj", 
                                   "cna1D_logFC", "cna1D_padj"),
                     col_types = c("text", "numeric", "numeric", 
                                   "numeric", "numeric")) %>%
  filter(!is.na(crz1D_logFC))


#+ load_af_data, echo=TRUE
af_ca_up.df = read_excel(eve.excelfile, sheet = 2, col_names = c("gene"), skip=1)
head(af_ca_up.df)
tail(af_ca_up.df)

af_ca_down.df = read_excel(eve.excelfile, sheet = 3, col_names = c("gene"), skip=1)
head(af_ca_down.df)
tail(af_ca_down.df)


af_crz_down.df = read_excel(eve.excelfile, sheet = 4, col_names = c("gene"), skip=1)
head(af_crz_down.df)
tail(af_crz_down.df)

af_crz_up.df = read_excel(eve.excelfile, sheet = 5, col_names = c("gene"), skip=1)
head(af_crz_up.df)
tail(af_crz_up.df)

#+ load_sc_data, echo=TRUE
sc_cna_ca.df = read_excel(eve.excelfile, sheet = 6, col_names = c("gene"), skip=2) %>%
  mutate(gene=str_extract(gene, pattern="\\w+")) %>% 
  filter(!is.na(gene))
head(sc_cna_ca.df)
tail(sc_cna_ca.df)

sc_cna_ca_na.df = read_excel(eve.excelfile, sheet = 7, col_names = c("gene"), skip=2)
head(sc_cna_ca_na.df)
tail(sc_cna_ca_na.df)

