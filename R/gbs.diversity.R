#' GBS Diversity
#' 
#' Calculate different measures of diveristy in a hap object
#' 
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' 
#' @param hap the hap object to manipulate
#' 
#' @keywords
#' 
#' @examples
#' 
#' @export

gbs.diversity <- function(hap){
  
  # TODO Calculate maf or error if absent
  
  tHap = hap
  a = substr(hap$alleles,1,1)
  b = substr(hape$alleles,3,3)
  
  # Recalculate allele counts
  tHap$alleleA = rowSums(tHap[,2:ncol(tHap)] == a, na.rm = T) 
  tHap$alleleB = rowSums(tHap[,2:ncol(tHap)] == b, na.rm = T) 
  tHap$missing = rowSums(tHap[,2:ncol(tHap)] == "N", na.rm = T) 
  tHap$het = rowSums(tHap[,2:ncol(tHap)] == "N", na.rm = T) 
  tHap$maf = apply(cbind(tHap$alleleA, tHap$alleleB), 1, min)/ apply(cbind(tHap$alleleA, tHap$alleleB, tHap$het), 1, sum)
  
  # (1 - p^2 - q^2 )/ # loci
  # Do na values need to be removed before calculating mean to get correct length? A vs B
  
  # A
  nei = mean(1-hap$maf^2-(1-hap$maf)^2,na.rm=T)
  
  # B
  nei = (1-hap$maf^2-(1-hap$maf)^2)
  nei = mean(nei[!is.na(nei)])
  

  # FST, possibly use ‘hierfstat’ package
}