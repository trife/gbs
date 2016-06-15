#' Read hap file(s) or data frame.
#'
#' Read a file or data frame to create a gbs object
#'
#' @author Trevor Rife, \email{trife@@ksu.edu}
#'
#' @param hap A data frame, file, or folder containing hap files.
#' @param delim The delimter used in the file(s).
#' @param data.col The column number for the first line or individual.
#'
#' @details
#' This function creates a gbs object from hap files or a data frame. The gbs object contains a data frame representing the header columns, a data frame representing the genotypic calls, an empty geno data frame, and an empty stats data frame. Both the geno and stats data frames are created with the \code{gbs.summary} function. If a folder is passed to \code{hap.read}, \code{hap.join} is called to merge all files within the folder and remove duplicate tags.
#'
#' @return A gbs object.
#'
#' @examples
#' data(wheat)
#' hap = hap.read(raw.hap)
#'
#' @export

hap.read <- function(hap.obj, delim="\t", data.col) {

  if(class(hap.obj)=="gbs") {
    cat(substitute(hap.obj), "is already of type gbs.")
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt))
    stop()
  }

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
      hap <- data.table::fread(input=hap.obj, sep=delim, check.names=FALSE, header=TRUE, data.table=F, strip.white=T)
    }
  } else {
    hap <- hap.obj
  }

  row.names(hap) <- 1:nrow(hap)
  header <- hap[,1:(data.col-1)]
  calls <- hap[,data.col:ncol(hap)]

  # Set class
  output <- list(header=header,calls=calls,geno=matrix(),stats=data.frame())
  class(output) <- "gbs"
  invisible(output)
}
