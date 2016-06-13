#' Allele match.
#'
#' Compares alleles across different individuals in a hap object.
#'
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#' @author Trevor Rife, \email{trife@@ksu.edu}
#'
#' @param hap A matrix or data frame consisting consisting of rows (markers) and columns (individuals).
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
  
  if(!"count"%in%result && !"percent"%in%result) {
    stop("Result type must be specified.")
  }

  allele.match <- hap

  cat('\n')
  message("Supplied dataset has ", nrow(hap), " SNPs and ", ncol(hap), " individuals.")
  message('Computing percent identity...')
  cat('\n')

  nS <- ncol(allele.match)
  id <- matrix(NA, nrow=nS, ncol=nS)

  pb <- txtProgressBar(min = 0, max = nrow(id), initial = 0, style = 3)

  for (i in 1:nrow(id)){
    id_pc <- rep(NA, length(i:nrow(id)))
    id_ct <- rep(NA, length(i:nrow(id)))
    line1 <- as.character(allele.match[,i])
    for (j in i:ncol(id)){
      line2 <- as.character(allele.match[,j])
      shared <- line1!="N" & line2!="N" & line1!="H" & line2!="H" # TODO redo for IUPAC
      common <- line1[shared] == line2[shared]
      if(any(result == "percent")) {
        id_pc[j-i+1] <- sum(common, na.rm = T)/sum(shared, na.rm = T)
      }

      if(any(result == "count")) {
        id_ct[j-i+1] <- sum(shared, na.rm = T)
      }
    }
    id[i:ncol(id),i] <- id_ct
    id[i,i:ncol(id)] <- id_pc

    setTxtProgressBar(pb,i)
  }

  rownames(id)=colnames(allele.match)
  colnames(id)=colnames(allele.match)

  if(graph) {
    par(mfrow=c(1,2))
    
    if(any(result == "percent")) {
      hist(id[upper.tri(id)]*100, xlab = paste('%', 'identity'), main = 'Distribution of % identity')
    }
    
    if(any(result == "count")) {
      hist(id[lower.tri(id)], xlab = '# of comparisons', main = 'Distribution of # comparisons')
    }
    
    par(mfrow=c(1,1))
  }

  return(id)
}
