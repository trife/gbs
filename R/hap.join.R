#' Hap Join
#' 
#' Combines multiple hap files into a single hap object
#' 
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' 
#' @param hap.dir directory where hap files are located
#' 
#' @keywords
#' 
#' @examples
#' 
#' @export

hap.join <- function(hap.dir,delim="\t"){
  joined = NULL
  hap = NULL
  fileList = list.files(path=hap.dir,full.names=TRUE)
  
  for(i in 1:length(fileList))
  {
    curFile = fileList[i]
    d = read.delim(file=curFile, header=TRUE, sep="\t", row.names=NULL,check.names=FALSE)
    cnHAP = colnames(d)[2:length(colnames(d))]
    d = d[,1:dim(d)[2]-1]
    colnames(d)=cnHAP
    d$QCcode = as.numeric(as.character(d$QCcode))
    print(paste(nrow(d)," in file ",i,sep=""))
    d$center=i
    joined = rbind(joined,d)
    
  }
  
  print(paste(nrow(joined)," markers total.",sep=""))
  
  #TODO print stats
  #TODO check data to verify
  #TODO print out stats about dfs
  
  invisible(joined)
}