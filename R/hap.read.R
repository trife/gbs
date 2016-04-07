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

hap.read <- function(file,delim="\t",data){
  #TODO if folder passed, use hap.join, else read in data

  hap = fread(file=file, sep=delim,check.names=FALSE,header=TRUE, data.table = F)

  #TODO check data


  invisible(hap)
}
