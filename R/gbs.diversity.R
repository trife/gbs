#' GBS diversity.
#'
#' Calculates different measures of diveristy in a gbs object.
#'
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#' @author Trevor Rife, \email{trife@@ksu.edu}
#'
#' @param hap The gbs object to manipulate.
#' @param clusters A data frame with two columns, where first column contain the names of individuals and second column contain the corresponding group (numeric) of population to which each individual belongs to.
#' @param maf.thresh
#' @param graph A logical option to graph the percent polymorphism.
#' @param het The symbol(s) used for heterozygous calls.
#' @param missing The symbol(s) used for missing calls.
#'
#' @details
#'
#' @return A list containing the Nei's genetic diversity and FST of the individuals
#'
#' @examples
#' data(wheat)
#' gbs.diversity(hap,graph=T)
#'
#' @export

gbs.diversity <- function(hap, clusters=NULL, maf.thresh=0.05, graph=F, het="H", missing="N") {

  if(class(hap)!="gbs") {
    stop("hap argument is incorrect type. Use hap.read to create gbs object.")
  }

  output = list()

  a <- substr(hap$header$alleles,1,1)
  b <- substr(hap$header$alleles,3,3)

  # Recalculate allele counts
  missing <- rowSums(hap$calls == missing, na.rm = T)
  het <- rowSums(hap$calls == het, na.rm = T)
  alleleA <- rowSums(hap$calls == a, na.rm = T)*2 + het
  alleleB <- rowSums(hap$calls == b, na.rm = T)*2 + het

  # Minor allele frequency
  maf <- apply(cbind(alleleA, alleleB), 1, min) / apply(cbind(alleleA, alleleB), 1, sum)

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
  nei <- mean(1-(maf^2)-(1-maf)^2, na.rm = T)
  cat("Nei's diversity index:", nei)
  output$nei = nei

  # TODO FST
  if(is.null(clusters)) {
    cat("Population groups not defined. Unable to calculate FST.")
  } else {
     tHapCalls <- t(hap$calls) # transformation required - ind in rows, markers in cols
     popGroups <- clusters[match(rownames(tHapCalls), clusters[,1]), ] # ordering the ind in cluster according to hap
     genidObject<-df2genind(tHapCalls, ncode = 1)

     cat("Computing overall population Fst...")
     popFst <- fstat(genidObject, pop = popGroups[,2], fstonly = T)
     output$popFst <- popFst
     cat(' Done')

     cat("computing Nei's (1973) pairwise Fst...")
     pairwiseFst <- pairwise.fst(genidObject, pop = popGroups[,2]) # according to Nei (1973)
     output$pairwiseFst <- pairwiseFst
     cat(' Done')
  }

  invisible(output)
}
