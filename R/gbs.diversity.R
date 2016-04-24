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

gbs.diversity <- function(hap, mafThresh=0){

  # TODO Calculate maf or error if absent

  tHap = hap
  a = substr(tHap$alleles,1,1)
  b = substr(tHap$alleles,3,3)

  # Recalculate allele counts
  tHap$missing = rowSums(tHap == "N", na.rm = T)
  tHap$het = rowSums(tHap == "H", na.rm = T)
  tHap$alleleA = rowSums(tHap == a, na.rm = T)*2 + tHap$het
  tHap$alleleB = rowSums(tHap == b, na.rm = T)*2 + tHap$het

  # minor allele frequency
  tHap$maf = apply(cbind(tHap$alleleA, tHap$alleleB), 1, min) / apply(cbind(tHap$alleleA, tHap$alleleB), 1, sum)

  # proportion of polymorphism
  par(mfrow=c(1,2))
  polyMarkers = sum(tHap$maf >= mafThresh)
  totalMarkers = nrow(tHap)
  cat('Number of polymorphic sites:', polyMarkers)
  cat('Proportion of polymorphic sites (P):', polyMarkers / totalMarkers)
  barplot(c(polyMarkers, totalMarkers-polyMarkers), xlim = c(0,5), col = c('lightgray','black'), legend.text = c('Polymorphic','Monomorphic'), ylab = "# markers")
  pie(c(polyMarkers, totalMarkers-polyMarkers), labels = c('Polymorphic markers','Monomorphic markers'))
  dev.off()

  # (1 - p^2 - q^2 )/ # loci
  # Do na values need to be removed before calculating mean to get correct length? A vs B

  # A
  nei = mean(1-hap$maf^2-(1-hap$maf)^2,na.rm=T)

  # B
  nei = (1-hap$maf^2-(1-hap$maf)^2)
  nei = mean(nei[!is.na(nei)])



  # FST, possibly use ‘hierfstat’ package
}
