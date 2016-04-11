install.packages("devtools")
install.packages("roxygen2")

library("devtools")
library(roxygen2)

document(pkg = ".")

install(".")
library("gbs")

help(allele.match)

install_github("trife/gbs")

# R CMD Rd2pdf gbs/