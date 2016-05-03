#' Bind multiple hap files together
#' 
#' Combines multiple hap files into a single hap object.
#' 
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' 
#' @param hap.dir The directory where the hap files are located.
#' @param delim The file delimiter for the hap files.
#' 
#' @seealso \code{\link{hap.read}}
#'
#' @export

hap.join <- function(hap.dir,delim="\t"){
  joined = NULL
  hap = NULL
  fileList = list.files(path=hap.dir,full.names=TRUE)
  
  for(i in 1:length(fileList)) {
    curFile = fileList[i]
    d = data.table::fread(input=curFile, header=TRUE, sep=delim, check.names=FALSE,data.table = F,strip.white=T)
    print(paste(nrow(d)," in file ",i,sep=""))
    d = cbind(dif=i,d)
    joined = rbind(joined,d)
  }
  
  print(paste(nrow(joined)," markers total.",sep=""))
  
  invisible(joined)
}