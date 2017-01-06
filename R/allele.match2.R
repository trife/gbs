#' Allele match.
#'
#' Compares alleles across different individuals in a hap object.
#'
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#' @author Trevor Rife, \email{trife@@ksu.edu}
#'
#' @param hap An object of type gbs.
#' @param result The values to be returned. Either the number of matched calls (\code{count}), \code{percent} identitiy, or both.
#' @param histgraph Logical value to graph results.
#'
#' @details
#' This function ignores missing data and heterozygous calls.
#'
#' @return A matrix with upper triangle equal to the percent identity of the individuals and lower triangle equal to the proportion of alleles that match between two individuals.
#'
#' @examples
#' data(wheat)
#' percId = allele.match(hap$calls,graph=T)
#'
#' @export

allele.match <- function(hap, result=c("count","percent"), graph=F) {

   if(is.null(dim(hap$calls))) {
      stop("Supplied dataset is empty!")
   }

  if(class(hap)!="gbs") {
    stop("Supplied dataset object is not of type gbs. Please use hap.read() to create gbs object first.")
  }

  if(!"count"%in%result && !"percent"%in%result) {
    stop("A valid result type must be specified.")
  }

  allele.match <- hap$calls

  cat('\n')
  message("Supplied dataset has ", nrow(allele.match), " SNPs and ", ncol(allele.match), " individuals.")
  message('Computing percent identity...')
  cat('\n')

  allele.match[allele.match=='H' | allele.match=='N']=NA
  nS <- ncol(allele.match)
  id <- matrix(NA, nrow=nS, ncol=nS)

  pb <- txtProgressBar(min = 0, max = nrow(id), initial = 0, style = 3)

  #percent id and number of comparisons computation for n-1 individuals
  for (i in 1 : (nrow(id)-1)){
     refInd <- allele.match[,i] #reference individual to be match against
     toBeMatched <- allele.match[,i:nS] #individuals to be match against refInd

        numComp <- colSums(!is.na(refInd == toBeMatched), na.rm = T) #total number of comparisons
        pcId <- colSums(refInd == toBeMatched, na.rm = T)/numComp #percent identity

        id[i:nS, i] = numComp
        id[i ,i:nS] = pcId

     setTxtProgressBar(pb,i+1)
  }

  #percent id for the last individual
  finalInd <- allele.match[,nS]
  numComp <- sum(!is.na(finalInd == finalInd), na.rm = T)
  pcId <- sum(finalInd == finalInd, na.rm = T) / numComp

  id[nS, nS] = pcId

  rownames(id)=colnames(allele.match)
  colnames(id)=colnames(allele.match)

  # if(graph) {
  #   par(mfrow=c(1,2))
  #
  #   if(any(result == "percent")) {
  #     hist(id[upper.tri(id)]*100, xlab = paste('%', 'identity'), main = 'Distribution of % identity')
  #   }
  #
  #   if(any(result == "count")) {
  #     hist(id[lower.tri(id)], xlab = '# of comparisons', main = 'Distribution of # comparisons')
  #   }
  #
  #   par(mfrow=c(1,1))
  # }

  return(id)
}
