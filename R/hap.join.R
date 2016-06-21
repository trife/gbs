#' Bind multiple hap files together
#'
#' Combines multiple hap files into a single hap object.
#'
#' @author Trevor Rife, \email{trife@@ksu.edu}
#'
#' @param hap.dir The directory where the hap files are located. Folder should only contain the hap files.
#' @param delim The file delimiter for the hap files.
#'
#' @details
#' For use with the old TASSEL 4 GBS pipeline that produced multiple files.
#'
#' @seealso \code{\link{hap.read}}
#'
#' @export

hap.join <- function(hap.dir,delim="\t"){
  joined <- NULL
  hap <- NULL
  fileList <- list.files(path=hap.dir,full.names=TRUE)

  for(i in 1:length(fileList)) {
    curFile <- fileList[i]
    d1 <- data.table::fread(input=fileList[1], header=TRUE, sep=delim, check.names=FALSE,data.table = F,strip.white=T)
    di <- data.table::fread(input=curFile, header=TRUE, sep=delim, check.names=FALSE,data.table = F,strip.white=T)
    print(paste(nrow(di)," rows in ", substitute(curFile),sep=""))
    di$center=i

    if(colnames(di)!=colnames(d1)) {
      stop("Column names in hap files do not match.")
    } else {
      joined <- rbind(joined,di)
    }
  }

  # Get rid of duplicates
  rs_pos <- cbind(joined$rs,joined$assembly)
  dup <- duplicated(rs_pos)
  noDup <- joined[!dup,]
  odr <- order(as.vector(noDup$rs), as.vector(noDup$assembly))
  joined <- noDup[odr,]
  joined$rs <- paste("rs", c(1:nrow(joined)), sep="")

  print(paste(nrow(joined)," markers total.",sep=""))

  invisible(joined)
}
