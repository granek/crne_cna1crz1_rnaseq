library(dplyr)
library(stringr)
library(tools)

af_file = "rewiring_analysis/fold_change_tables/soriani_2008.csv"
af.df = read.csv(af_file, stringsAsFactors = FALSE, sep=",", skip=1, 
                 strip.white=TRUE, col.names = c("desc","log2fc")) %>%
  mutate(log2fc=as.numeric(log2fc)) %>%
  filter(!is.na(log2fc)) %>%
  mutate(gene_name = str_sub(desc,1,10)) %>%
  select(gene_name, log2fc)

clean_af_file = paste0(file_path_sans_ext(af_file),"_clean.csv")
write.csv(af.df,file=clean_af_file,row.names = FALSE,quote=FALSE)
