#' Read hap file or data frame
#'
#' Read in a file or data frame to create a hap object
#'
#' @author Trevor Rife, \email{trife@@ksu.edu}
#'
#' @param hap the hap object, file, or folder
#' @param delim the delimter used in the file(s)
#' @param data column number for the first individual. Removes unnecessary columns from the data frame if specified
#' @param genotypes list of calls used in the hap object
#'
#' @keywords
#'
#' @examples
#'
#' @export

hap.read <- function(hap.obj, delim="\t", data, genotypes=c(NA,"A","T","C","G","H","N")){
  if(!file.exists(hap.obj)) {
    stop("File or folder does not exist.")
  }

  if(file.info(hap.obj)$isdir) {
    hap = hap.join(hap.obj,delim)
  }
  
  if(!file.info(hap.obj)$isdir) {
    hap = data.table::fread(input=hap.obj, sep=delim,check.names=FALSE,header=TRUE, data.table = F,strip.white=T)
    hap = cbind(dif=1,hap)
  }

  if(!"rs"%in%colnames(hap)) {
    stop("Tag column (rs) missing from hap object")
  }
  
  if(!"alleles"%in%colnames(hap)) {
    stop("Alleles column (alleles) missing from hap object")
  }
  
  if(!"pos"%in%colnames(hap)) {
    stop("Marker position column (pos) missing from hap object")
  }
  
  if(!"chrom"%in%colnames(hap)) {
    stop("Chromosome column (chrom) missing from hap object")
  }
  
  if(!missing(data)) {
    hap = cbind(dif=hap$dif,rs=hap$rs,alleles=hap$alleles,pos=hap$pos,chrom=hap$chrom,hap[,(data+1):ncol(hap)])
    
    # Simple data check for non-genotypes in first 10 roww
    if(!all(apply(hap[1:10,6:ncol(hap)],MARGIN=2,function(x) x%in%genotypes))) {
      stop("Non genotypes detected in hap object. Edit the genotypes parameter or check your hap object.")
    }
  }
  
  # Simple data check for non-genotypes in first 10 rows
  if(!all(apply(hap[1:10,6:ncol(hap)],MARGIN=2,function(x) x%in%genotypes))) {
    stop("Non genotypes detected in hap object. Edit the genotypes parameter or check your hap object.")
  }
  
  invisible(hap)
}