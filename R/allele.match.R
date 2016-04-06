#' Allele match
#' 
#' This function compares alleles across lines in a hap object
#' 
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#' 
#' @param hap a hap object
#' 
#' @keywords alleles, hap
#' 
#' @examples
#' 
#' @export

allele.match <- function(hap){
  hap.obj = hap
  allele.match = hap.obj[,-c(1:11)]
  allele.match[allele.match=="H" | allele.match=="N"]=NA
  dim(allele.match)
  allele.match[1:5,1:20]
  
  nS = ncol(allele.match)
  id = matrix(NA, nrow=nS, ncol=nS)
  
  ## New
  for (i in 1:nrow(id)){
    id_pc <- rep(NA, length(i:nrow(id)))
    id_ct <- rep(NA, length(i:nrow(id)))
    line1 = as.character(allele.match[,i])
    for (j in i:ncol(id)){
      line2 = as.character(allele.match[,j])
      shared = line1!="N" & line2!="N" & line1!="H" & line2!="H"
      common = line1[shared] == line2[shared]
      id_pc[j-i+1] <- sum(common, na.rm = T)/sum(shared, na.rm = T)
      id_ct[j-i+1] <- sum(shared, na.rm = T)
    }
    id[i:ncol(id),i] <- id_ct
    id[i,i:ncol(id)] <- id_pc
  }
  
  
  rm(common, id_ct, id_pc, line1, line2, shared)
  
  rownames(id)=colnames(allele.match)
  colnames(id)=colnames(allele.match)
  id[1:10,1:10]
  
  #TODO change output types
}