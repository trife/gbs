#' GBS Graphs
#' 
#' Calculates and graphs relevant stats for a hap object
#' 
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' 
#' @param hap the hap object to be processed
#' 
#' @keywords
#' 
#' @examples
#' 
#' @export


gbs.graph <- function(hap) {
  ## Minor allele frequency
  MAF = apply(cbind(hap$alleleA, hap$alleleB), 1, min)/ apply(cbind(hap$alleleA, hap$alleleB, hap$het), 1, sum)
  hist(MAF, xlab="minor allele freq", ylab="# snps")
  
  ## Percent heterozygous
  percentHET = hap$het / apply(cbind(hap$alleleA, hap$alleleB, hap$het), 1, sum)
  hist(percentHET)
  hist(hap$het)
  het = hap$het / (hap$alleleA + hap$alleleB + hap$het)
  hist(het, xlab="minor allele freq", ylab="# snps")
  
  ## Add to hap
  hap = cbind(hap[,c(1:8)], MAF, percentHET, hap[,c(9:ncol(hap))])
  
  ## Percent present
  hap$present = as.numeric(as.character(hap$present))
  sum(hap$present>=0.8)
  hist(hap$present,56)
  dim(hap)
  
  ## SNPs per line
  hapNA = hap=="N"
  hapNA[1:3,]
  lineData = apply(hap!="N", 2, sum)[-c(1:13)]
  hist(lineData, main=project, xlab="# SNPs", ylab="# lines",18)
}