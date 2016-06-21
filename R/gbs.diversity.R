#' GBS diversity.
#'
#' Calculates different measures of diveristy in a gbs object.
#'
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#' @author Trevor Rife, \email{trife@@ksu.edu}
#'
#' @param hap The gbs object to manipulate.
#' @param popGroups A data frame with two columns, where first column contain the names of individuals and second column contain the corresponding group (numeric) of population to which each individual belongs to.
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

gbs.diversity <- function(hap, popGroups=NULL, maf.thresh=0.05, graph=F, het="H", missing="N") {

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
    cat('Number of polymorphic sites with MAF threshold:', polyMarkers,"\n")
    cat('Proportion of polymorphic sites (P):', polyMarkers / totalMarkers,"\n")
    barplot(c(polyMarkers, totalMarkers-polyMarkers), xlim = c(0,5), col = c('lightgray','black'), legend.text = c('Polymorphic','Monomorphic'), ylab = "# markers")
    pie(c(polyMarkers, totalMarkers-polyMarkers), labels = c('Polymorphic \nmarkers','Monomorphic \nmarkers'))
  }

  # Nei's diversity index
  nei <- mean(1-(maf^2)-(1-maf)^2, na.rm = T)
  cat("Computing Nei's diversity index... Done")
  output$nei = nei

  # Fst (Population Fixation index)
  if(is.null(popGroups)) {
    cat("Population groups not defined. Unable to calculate Fst.")
  } else {
     # Numericalizing data
     cat('\n')
     cat('Numericalizing genotypes...')

     hapCalls=as.matrix(hap$calls)
     hapCalls[hapCalls=='N']=NA
     lst=apply(hapCalls, 1, function(x) {
        lev <- unique(x)
        lev <- lev[lev!='H' & !is.na(lev)]
        x[x=='H']=12
        x[x==lev[1]]=11
        x[x==lev[2]]=22
        return(x)
      }
     )
     cat(' Done')

     popGroups <- popGroups[match(rownames(lst), popGroups[,1]), ] # ordering the ind in cluster according to hap
     genidObject<-df2genind(lst, ncode = 2, NA.char = NA, sep = '/')


     cat('\n')
     cat("Computing overall population Fst...")
     popFst <- fstat(genidObject, pop = popGroups[,2], fstonly = T)
     output$popFst <- popFst
     cat(' Done')

     cat('\n')
     cat("computing Nei's (1973) pairwise Fst...")
     pairwiseFstNei1973 <- pairwise.fst(genidObject, pop = popGroups[,2], res.type = "matrix") # according to Nei (1973)
     output$pairwiseFstNei1973 <- pairwiseFstNei1973
     cat(' Done')

     dat=as.data.frame(cbind(popGroups[,2],lst))
     for (i in 1:ncol(dat)) {
        dat[,i]=as.character(dat[,i])
        dat[,i]=as.numeric(dat[,i])
     }

     cat('\n')
     cat("computing Nei's (1987) pairwise Fst...")
     pairwiseFstNei1987 <- pairwise.neifst(dat = dat)
     output$pairwiseFstNei1987 <- pairwiseFstNei1987
     cat(' Done')

     cat('\n')
     cat("computing Weir and Cockerham (1984) pairwise Fst...")
     pairwiseFstWC <- pairwise.WCfst(dat = dat)
     output$pairwiseFstWC <- pairwiseFstWC
     cat(' Done')

   }

  invisible(output)
}
