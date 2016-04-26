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

hap.read <- function(hap, delim="\t", data, genotypes=c(NA,"A","T","C","G","H","N")){
  if(!file.exists(file)) {
    stop("File or folder does not exist.")
  }

  if(file.info(file)$isdir) {
    hap = hap.join(file,delim)
  }
  
  if(!file.info(file)$isdir) {
    hap = data.table::fread(input=file, sep=delim,check.names=FALSE,header=TRUE, data.table = F,strip.white=T)
    hap = cbind(dif=1,hap)
  }
  
  if(is.data.frame(hap)) {
    # Simple data check for non-genotypes in first 10 rows
    if(!all(apply(hap[1:10,5:ncol(hap)],MARGIN=2,function(x) x%in%genotypes))) {
      stop("Non genotypes detected in hap object. Edit the genotypes parameter or check your hap object.")
    }
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
    hap = cbind(hap[,"rs"],hap[,"alleles"],hap[,"pos"],hap[,"chrom"],hap[,data:ncol(hap)])
  }
  
  invisible(hap)
}