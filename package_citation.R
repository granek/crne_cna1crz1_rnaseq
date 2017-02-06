bib_file = "r_packages.bib"
file.remove(bib_file)

package_vec = c("gplots","RColorBrewer","optparse")

for (pack in package_vec) {
  print(pack)
  library(pack,character.only = TRUE)
  write(toBibtex(citation(pack)),file=bib_file,append=TRUE)
}
