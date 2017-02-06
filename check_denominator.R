library(dplyr)
library(reshape2)
library(tidyr)

setwd("~/Documents/BioinfCollabs/HeitmanLab/calcineurin_reg/make_cna/results")

count_file = "cna1_crz1_allcounts.csv"
counts = read.csv(count_file,row.names=1)

##=======================================================================
##=======================================================================
results_file = "condition_KO_cna1_vs_WT_37C_results_5fdr_2fc.csv"
grep_strains = "cna|WT"
grep_temp = "37C"

results = read.csv(results_file,row.names=1) %>% add_rownames %>% arrange(log2FoldChange) %>% select(rowname,log2FoldChange)
extremes = rbind(head(results,2),tail(results,2))

for (gene in extremes$rowname) {
  gene_counts = melt(counts[gene,]) %>% filter(grepl(grep_temp,variable)) %>% filter(grepl(grep_strains,variable)) %>% arrange(value)
  print(gene)
  print(gene_counts)
}
print(extremes)
##=======================================================================
##=======================================================================
results_file = "condition_KO_crz1_vs_WT_37C_results_5fdr_2fc.csv"
grep_strains = "crz|WT"
grep_temp = "37C"

results = read.csv(results_file,row.names=1) %>% add_rownames %>% arrange(log2FoldChange) %>% select(rowname,log2FoldChange)
extremes = rbind(head(results,2),tail(results,2))

for (gene in extremes$rowname) {
  gene_counts = melt(counts[gene,]) %>% filter(grepl(grep_temp,variable)) %>% filter(grepl(grep_strains,variable)) %>% arrange(value)
  print(gene)
  print(gene_counts)
}
print(extremes)
##=======================================================================
##=======================================================================
##=======================================================================
# OLD Stuff
## separate(key, into = c("location", "time"), sep = "\\.") 

