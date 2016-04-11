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
#' @keywords hap
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
    hap = read.delim(file=file, sep=delim,check.names=FALSE,header=TRUE)
  }
  
  #TODO check data
  
  
  invisible(hap)
}