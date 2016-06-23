#' Bind multiple hap files together
#'
#' Combines multiple hap files into a single hap object.
#'
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' @author Narinder Singh, \email{nss470@@ksu.edu}
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

hap.join <- function(hap.dir,delim="\t",data.col){

   cat('\n')

   if(missing(data.col)) {
      stop("Specify which column contains the first individual.")
   }

  joined <- NULL
  fileList <- list.files(path=hap.dir,full.names=TRUE)

  for(i in 1:length(fileList)) {
    curFile <- fileList[i]
    hapDirPathLength = length(strsplit(curFile, "/")[[1]])
    hapFileName = strsplit(curFile, "/")[[1]][hapDirPathLength]

    cat('File "', substitute(hapFileName), '" has... ', sep = "")
    di <- data.table::fread(input=curFile, header=TRUE, sep=delim, check.names=FALSE,data.table = F,strip.white=T)

    if(i==1) d1=di

    cat(nrow(di), ' rows', sep = ""); cat('\n')
    di$center=i

    if(!all(colnames(di)==colnames(d1))) {
      stop('Column names in "', hapFileName, '" do not match the column names in other files.')
    } else {
      joined <- rbind(joined,di)
    }
  }

  # Get rid of duplicates
  cat('Removing duplicate tags... ')
  rs_pos <- cbind(joined$rs,joined$assembly)
  dup <- duplicated(rs_pos)
  noDup <- joined[!dup,]
  odr <- order(as.vector(noDup$rs), as.vector(noDup$assembly))
  joined <- noDup[odr,]
  cat('Done'); cat('\n')

  joined <- hap.read(joined, data.col = data.col)

  cat('\n'); cat('\n')
  invisible(joined)
}
