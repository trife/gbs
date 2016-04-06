#' GBS Diversity
#' 
#' Calculate different measures of diveristy in a hap object
#' 
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' 
#' @param hap the hap object to manipulate
#' @param delim the delimiter for the hap object (if file)
#' 
#' @keywords
#' 
#' @examples
#' 
#' @export

gbs.diversity <- function(hap,delim="\t"){
  require(ggplot2)
  require(rrBLUP)
  require(BLR)
  require(multicore)
  require(parallel)
  require(ape)
  
  hap.obj = read.table(file=paste("R/",project,"/hap/",project,"_srpn_hap.txt",sep=""), header=TRUE,sep="\t",check.names=FALSE)
  hap.obj[1:20,1:50]
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
  
  ## Write and read table (converts to numeric)
  write.table(hap01, file=paste("R/",project,"/hap/",project,"_srpn_hap01.txt",sep=""), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
  hap01 = read.table(file=paste("R/rpn/hap/rpn_srpn_hap01.txt",sep=""), header=TRUE,check.names=FALSE,sep="\t")
  
  ## Calculate A matrix using rrBLUP
  geno = as.matrix(hap01[,12:ncol(hap01)])
  geno = t(geno)
  geno[1:20,1:20]
  dim(geno)
  write.table(geno,file="R/rpn/hap/rpn_srpn_geno.txt",sep="\t",col.names=TRUE,row.names=TRUE,quote=F)
  geno = read.table(file="R/rpn/hap/rpn_srpn_geno.txt",sep="\t")
}