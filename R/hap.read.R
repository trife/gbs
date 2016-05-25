#' Read hap file(s) or hap object.
#'
#' Read a file or data frame to create a hap object
#'
#' @author Trevor Rife, \email{trife@@ksu.edu}
#'
#' @param hap The hap object, file, or folder containing hap files.
#' @param delim The delimter used in the file(s).
#' @param data.col The column number for the first line or individual.
#'
#' @details
#' If a folder is passed to \code{hap.read}, \code{hap.join} is called to merge all files within the folder and remove duplicates.
#'
#' @return A verified hap data frame.
#'
#' @examples
#' data(wheat)
#' hap = hap.read(wheat)
#'
#' @export

hap.read <- function(hap.obj, delim="\t", data.col) {
  
  if(missing(data.col)) {
    stop("Specify which column contains the first individual.")
  }
  
  if(class(hap.obj) != "data.frame") {
    if(!file.exists(hap.obj)) {
      stop("File or folder does not exist.")
    }
    
    if(file.info(hap.obj)$isdir) {
      hap <- hap.join(hap.obj, delim)
    }
    
    if(!file.info(hap.obj)$isdir) {
      hap <- data.table::fread(input=hap, sep=delim, check.names=FALSE, header=TRUE, data.table = F, strip.white=T)
    }
  } else {
    hap = hap.obj
  }

  row.names(hap) <- 1:nrow(hap)
  header <- hap[,1:(data.col-1)]
  calls <- hap[,data.col:ncol(hap)]
  
  # Set class
  output <- list(header=header,calls=calls,geno=matrix(),stats=data.frame())
  class(output) <- "gbs"
  invisible(output)
}