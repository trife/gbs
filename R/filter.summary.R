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

filter.summary <- function(hap,project="gbs",output=c("hap","geno")){
  rs_pos_alleles = hap[,c(1,6,2)]
  dim(rs_pos_alleles)
  rs_pos = hap[,c(1,6)]
  
  dup = duplicated(rs_pos)

  noDup = hap[!dup,]
  as.vector(noDup$rs)
  odr = order(as.vector(noDup$rs), as.vector(noDup$assembly))
  
  hap = noDup[odr,]
 
  hap$pos = c(1:nrow(hap))
  rsOrig = hap$rs
  
  hap$rs = paste(project, c(1:nrow(hap)), sep="")
  hap = cbind(rsOrig,hap)
  
  ## Check blank wells for data
  if(any(grepl("BLANK",colnames(hap)))) {
    missing.blank = hap[,grepl("BLANK",colnames(hap))]=="N"
    blank = as.matrix(apply(!missing.blank, 2, sum))
    nrow(hap)
    blank
  }
  
  ## Rename some columns
  colnames(hap)[colnames(hap)=="assembly"] = "snp_pos"
  colnames(hap)[colnames(hap)=="protLSID"] = "alleleA"
  colnames(hap)[colnames(hap)=="assayLSID"] = "alleleB"
  colnames(hap)[colnames(hap)=="panelLSID"] = "het"
  colnames(hap)[colnames(hap)=="QCcode"] = "present"
  colnames(hap)[colnames(hap)=="center"] = "dif"
  
  ## Remove some columns
  hap = hap[,colnames(hap)!="chrom"]
  hap = hap[,colnames(hap)!="pos"]
  hap = hap[,colnames(hap)!="strand"]
  hap = hap[,!grepl("blank",colnames(hap), ignore.case=TRUE)]
  colnames(hap)
  
  #TODO output hap or geno
  if(any(output=="hap")) {
    invisible(hap)
  }
  if(any(output=="geno")) {
    hap01 = hap
    
    ## Define allele a and allele b
    a = substring(hap$alleles,1,1)
    a[hap$alleleA<hap$alleleB] = substring(hap$alleles,3,3)[hap$alleleA<hap$alleleB]
    b = substring(hap$alleles,3,3)
    b[hap$alleleA<hap$alleleB] = substring(hap$alleles,1,1)[hap$alleleA<hap$alleleB]
    sum(a == b)
    
    ## Turn allele a and allele b into -1 and 1.  Het into 0
    hap01[hap == a] = "-1"
    hap01[hap == b] = "1"
    hap01[hap == "H"] = "0"
    
    ## TODO convert to numeric
    
    
    ## Calculate A matrix using rrBLUP
    geno = as.matrix(hap01[,12:ncol(hap01)])
    geno = t(geno)
    invisible(geno)
  }
}