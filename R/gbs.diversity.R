#' Determine GBS diversity.
#'
#' Calculates different measures of diveristy in a hap object.
#'
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#' @author Trevor Rife, \email{trife@@ksu.edu}
#'
#' @param hap The hap object to manipulate.
#' @param clusters A data frame with row names equal to individuals and a single column dictating which group each individual belongs to.
#' @param maf.thresh
#' @param graph Logical option to graph percent polymorphism.
#'
#' @details
#'
#' @return A list containing the Nei's genetic diversity and FST of the individuals
#'
#' @examples
#' data(wheat)
#' hap = hap.read(wheat)
#' hap = filter.summary(hap,project="wheat",output="hap")$hap
#' gbs.diversity(hap,data.col=14,graph=T)
#'
#' @export

gbs.diversity <- function(hap, clusters=NULL, maf.thresh=0.05, graph=F){
  
  if(class(hap)!="gbs") {
    stop("hap argument is incorrect type. Use hap.read to create gbs object.")
  }
  
  output = list()
  
  a = substr(hap$header$alleles,1,1)
  b = substr(hap$header$alleles,3,3)
  
  # Recalculate allele counts
  missing = rowSums(hap$calls == "N", na.rm = T) # TODO change for IUPAC
  het = rowSums(hap$calls == "H", na.rm = T) # TODO change for IUPAC
  alleleA = rowSums(hap$calls == a, na.rm = T)*2 + het
  alleleB = rowSums(hap$calls == b, na.rm = T)*2 + het
  
  # Minor allele frequency
  maf = apply(cbind(alleleA, alleleB), 1, min) / apply(cbind(alleleA, alleleB), 1, sum)
  
  # Graph proportion of polymorphism
  if(graph) {
    par(mfrow=c(1,2))
    polyMarkers = sum(maf >= maf.thresh)
    totalMarkers = nrow(hap$calls)
    cat('Number of polymorphic sites:', polyMarkers,"\n")
    cat('Proportion of polymorphic sites (P):', polyMarkers / totalMarkers,"\n")
    barplot(c(polyMarkers, totalMarkers-polyMarkers), xlim = c(0,5), col = c('lightgray','black'), legend.text = c('Polymorphic','Monomorphic'), ylab = "# markers")
    pie(c(polyMarkers, totalMarkers-polyMarkers), labels = c('Polymorphic \nmarkers','Monomorphic \nmarkers'))
  }
  
  # Nei's diversity index
  nei = mean(1-(maf^2)-(1-maf)^2, na.rm = T)
  cat("Nei's diversity index:", nei,"\n")
  output$nei = nei
  
  # TODO FST
  if(is.null(clusters)) {
    cat("Clusters not specified; unable to calculate FST.")
  } else {
    output$fst = NULL
  }
  
  invisible(output)
}
