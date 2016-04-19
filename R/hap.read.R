#' Read hap file
#'
#' Read in a file to create a hap object
#'
#' @author Trevor Rife, \email{trife@@ksu.edu}
#'
#' @param file the hap file
#' @param delim the delimter used in the file
#' @param data column number where genotypes begin
#'
#' @keywords
#'
#' @examples
#'
#' @export

hap.read <- function(file,delim="\t"){
  if(!file.exists(file)) {
    stop("File or folder does not exist.")
  }

  # If folder passed, use hap.join, else read data
  if(file.info(file)$isdir) {
    hap = hap.join(file,delim)
  }

  if(!file.info(file)$isdir) {
    hap = data.table::fread(input=file, sep=delim,check.names=FALSE,header=TRUE, data.table = F,strip.white=T)
    hap = cbind(dif=1,hap)
  }

  if(!"rs"%in%colnames(hap)) {
    stop("Tag column missing (rs)")
  }
  
  if(!"alleles"%in%colnames(hap)) {
    stop("Alleles column (alleles) missing")
  }
  
  if(!"pos"%in%colnames(hap)) {
    stop("Marker position column (pos) missing")
  }
  
  # TODO add in break here based on structure
  
  invisible(hap)
}