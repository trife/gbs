#' GBS Diversity
#'
#' Calculate different measures of diveristy in a hap object
#'
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#' @author Trevor Rife, \email{trife@@ksu.edu}
#'
#' @param hap the hap object to manipulate
#' @param maf.thresh threshold for polymorphism
#' @param graph option to graph percent polymorphism
#'
#' @keywords
#'
#' @examples
#'
#' @export

gbs.diversity <- function(hap, maf.thresh=0, graph=F){

  output = list()
  
  tHap = hap
  
  if(!"alleles"%in%colnames(hap)) {
    stop("Alleles column (alleles) missing from hap object")
  }
  
  if(!"alleleA"%in%colnames(hap)) {
    stop("B allele column missing from hap object. Run filter.summary.")
  }
  
  if(!"alleleB"%in%colnames(hap)) {
    stop("A allele column missing from hap object. Run filter.summary.")
  }
  
  if(!"het"%in%colnames(hap)) {
    stop("het column missing from hap object. Run filter.summary.")
  }
  
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
  if(graph) {
    par(mfrow=c(1,2))
    polyMarkers = sum(tHap$maf >= maf.thresh)
    totalMarkers = nrow(tHap)
    cat('Number of polymorphic sites:', polyMarkers,"\n")
    cat('Proportion of polymorphic sites (P):', polyMarkers / totalMarkers,"\n")
    barplot(c(polyMarkers, totalMarkers-polyMarkers), xlim = c(0,5), col = c('lightgray','black'), legend.text = c('Polymorphic','Monomorphic'), ylab = "# markers")
    pie(c(polyMarkers, totalMarkers-polyMarkers), labels = c('Polymorphic \nmarkers','Monomorphic \nmarkers'))
  }

  # Nei's diversity index
  nei = mean(1-(tHap$maf^2)-(1-tHap$maf)^2, na.rm = T)
  cat("Nei's diversity index:", nei,"\n")
  output$nei = nei

  # FST, possibly use ‘hierfstat’ package
  
  invisible(output)
}
