#' Filter Summary
#' 
#' Processes hap object, calculates summary statistics, and graphs relevant stats
#' 
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' @author Jesse Poland, \email{jpoland@@ksu.edu}
#' 
#' @param hap the hap object to be processed
#' 
#' @keywords
#' 
#' @examples
#' 
#' @export

filter.summary <- function(hap,project="gbs",data.col=4,output=c("hap","geno")){
  #TODO assume we have chromosome
  #TODO get rid of the other columns that default with the file if they exist
  
  # Get rid of duplicates
  rs_pos = hap[,c(2,4)]
  dup = duplicated(rs_pos)
  noDup = hap[!dup,]
  
  # Reorder
  odr = order(as.vector(noDup$rs), as.vector(noDup$pos))
  hap = noDup[odr,]
  
  # Add new columns
  rsOrig = hap$rs
  hap$rs = paste(project, c(1:nrow(hap)), sep="")
  hap = cbind(rsOrig,hap)
  
  # Check blank wells for data
  if(any(grepl("blank",colnames(hap),ignore.case=TRUE))) {
    missing.blank = hap[,grepl("BLANK",colnames(hap),ignore.case=TRUE)]=="N"
    blank = as.matrix(apply(!missing.blank, 2, sum))
    nrow(hap)
    blank
    
    # Remove blank columns
    hap = hap[,!grepl("blank",colnames(hap), ignore.case=TRUE)]
  }
  
  # Calculate allele stats TODO this feels inefficient
  a = substring(hap$alleles,1,1)
  b = substring(hap$alleles,3,3)
  
  counts = apply(hap[,6:ncol(hap)],MARGIN=1,function(x) table(t(x)))
  
  alleleA = lapply(1:length(counts),function(x) counts[[x]][names(counts[[x]])==a[x]])
  alleleB = lapply(1:length(counts),function(x) counts[[x]][names(counts[[x]])==b[x]])
  het = lapply(1:length(counts),function(x) counts[[x]][names(counts[[x]])=="H"])
  missing = lapply(1:length(counts),function(x) counts[[x]][names(counts[[x]])=="N"])
  
  alleleA[!(sapply(alleleA, length))] = 0
  alleleB[!(sapply(alleleB, length))] = 0
  het[!(sapply(het, length))] = 0
  missing[!(sapply(missing, length))] = 0
  
  alleleA = unlist(alleleA)
  alleleB = unlist(alleleB)
  het = unlist(het)
  missing = unlist(missing)
  
  present = (alleleA + alleleB) / (alleleA+alleleB+missing+het)
  
  # Bind columns onto hap
  hap = cbind(hap[,1:5],alleleA,alleleB,het,present,hap[,6:ncol(hap)])
  
  # Calculate and add summary stats
  maf = apply(cbind(hap$alleleA, hap$alleleB), 1, min) / apply(cbind(hap$alleleA, hap$alleleB, hap$het), 1, sum)
  phet = hap$het / apply(cbind(hap$alleleA, hap$alleleB, hap$het), 1, sum)
  hap$present = as.numeric(as.character(hap$present))
  hap = cbind(hap[,c(1:8)], maf, phet, hap[,c(9:ncol(hap))])

  # Output to list
  hapReturn = list()

  if(any(output=="geno")) {
    hap01 = hap
    hap01[,12:ncol(hap01)]=NA
    
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
    
    ## Calculate A matrix using rrBLUP
    geno = as.matrix(hap01[,12:ncol(hap01)])
    class(geno) = "numeric"
    geno = t(geno)
    hapReturn$geno = geno
  }
  if(any(output=="hap")) {
    hapReturn$hap = hap
  }
  
  hapReturn
}

# Separate function to just convert to geno
filter.summary.geno <- function(hap,data.col=12) {
  hapReturn = list()
  
  hap01 = hap
  hap01[,data.col:ncol(hap01)]=NA
  
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
  
  ## Calculate A matrix using rrBLUP
  geno = as.matrix(hap01[,data.col:ncol(hap01)])
  class(geno) = "numeric"
  geno = t(geno)
  hapReturn$geno = geno
}
