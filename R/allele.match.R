#' Allele match
#'
#' This function compares alleles across lines in a hap object
#'
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#' @author Trevor Rife, \email{trife@@ksu.edu}
#'
#' @param hap a matrix or data frmae consisting consisting of rows (markers) and columns (individuals)
#' @param result return the number of matched calls, percent identitiy, or both
#' @param calls optional vector to include non-standard genotype calls
#'
#' @keywords alleles, hap
#'
#' @examples
#'
#' @export

allele.match <- function(hap,result=c("count","percent"),calls=NULL,histgraph=F){
  if(!"count"%in%result && !"percent"%in%result) {
    stop("Result type must be specified.")
  }

  allele.match <- hap
  uniqCalls <- unique(c(hap))

  # Check to ensure input matrix has correct format
  genotypes <- c(NA,"A","C","G","T","H","N","a","c","g","t","h","n")

  if(!missing(genotypes)) {
    genotypes <- c(genotypes, calls)
  }

  offCalls <- paste(shQuote(uniqCalls[!uniqCalls%in%genotypes]), collapse=", ")

  if(!all(apply(allele.match,MARGIN=2,function(x) x%in%genotypes))) {
    stop("Following non-genotypes calls detected in the input matrix: ", offCalls, ". Edit the calls parameter to allow these.")
  }

  cat('\n')
  message("Supplied dataset has ", nrow(allele.match), " SNPs and ", ncol(allele.match), " individuals.")
  message('Computing percent identity...')
  cat('\n')

  nS <- ncol(allele.match)
  id <- matrix(NA, nrow=nS, ncol=nS)

  pb = txtProgressBar(min = 0, max = nrow(id), initial = 0, style = 3)

  for (i in 1:nrow(id)){
    id_pc <- rep(NA, length(i:nrow(id)))
    id_ct <- rep(NA, length(i:nrow(id)))
    line1 = as.character(allele.match[,i])
    for (j in i:ncol(id)){
      line2 = as.character(allele.match[,j])
      shared = line1!="N" & line2!="N" & line1!="H" & line2!="H"
      common = line1[shared] == line2[shared]
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

  if(histgraph) {
     pdf(file = "identityHist.pdf", width = 11, height = 8.5)
     hist(id[upper.tri(id)]*100, xlab = paste('%', 'identity'), main = 'Distribution of % identity')
     hist(id[lower.tri(id)], xlab = '# of comparisons', main = 'Distribution of # comparisons')
     dev.off()
     cat('\n')
     message('Histograms are saved as "identityHist.pdf" in your working directory.')
     cat('\n')
  }

  save(id, file = 'idMatrix.RData')

  message('DONE computing successfully. Computed id matrix is saved as "idMatrix.RData" in your working directory for future usage and can be read in using the following command:')
  cat('\n')
  cat("\tload('idMatrix.RData')")
  return(id)
}
