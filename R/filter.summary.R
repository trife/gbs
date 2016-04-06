#' Filter Summary
#' 
#' Processes hap object, calculates summary statistics, and graphs relevant stats
#' 
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' 
#' @param hap the hap object to be processed
#' 
#' @keywords
#' 
#' @examples
#' 
#' @export

filter.summary <- function(hap){
  rs_pos_alleles = joined[,c(1,6,2)]
  dim(rs_pos_alleles)
  rs_pos_alleles[1:10,]
  rs_pos = joined[,c(1,6)]
  
  dup = duplicated(rs_pos)
  sum(dup)
  
  noDup = joined[!dup,]
  nrow(noDup)
  as.vector(noDup$rs)
  rank(noDup$rs)
  odr = order(as.vector(noDup$rs), as.vector(noDup$assembly))
  
  hap = noDup[odr,]
  hap[1:5,]
  dim(hap)
  unique(hap$rs)
  nrow(noDup)
  
  hap$pos = c(1:nrow(hap))
  rsOrig = hap$rs
  
  hap$rs = paste(project, c(1:nrow(hap)), sep="")
  hap = cbind(rsOrig,hap)
  
  ## Check blank wells for data
  missing.blank = hap[,grepl("BLANK",colnames(hap))]=="N"
  blank = as.matrix(apply(!missing.blank, 2, sum))
  nrow(hap)
  blank
  
  ## Rename some columns
  colnames(hap)[colnames(hap)=="assembly"] = "snp_pos"
  colnames(hap)[colnames(hap)=="protLSID"] = "alleleA"
  colnames(hap)[colnames(hap)=="assayLSID"] = "alleleB"
  colnames(hap)[colnames(hap)=="panelLSID"] = "het"
  colnames(hap)[colnames(hap)=="QCcode"] = "present"
  colnames(hap)[colnames(hap)=="center"] = "dif"
  
  colnames(hap)[1:20]
  colnames(hap)[7]
  hap[1:3,1:50]
  
  ## Remove some columns
  hap = hap[,colnames(hap)!="chrom"]
  hap = hap[,colnames(hap)!="pos"]
  hap = hap[,colnames(hap)!="strand"]
  hap = hap[,!grepl("blank",colnames(hap), ignore.case=TRUE)]
  #hap = hap[,!colSums(hap=="N")>44000]
  colnames(hap)
}