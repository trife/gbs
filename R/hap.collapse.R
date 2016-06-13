#' Merge duplicate individuals.
#' 
#' Combines multiple individuals that have the same name (e.g. ind.1,ind.2,...,ind.n) in a gbs object.
#' 
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' 
#' @param hap A gbs object consisting only of individuals (columns) and markers (rows).
#' @param names A list of the base names of the individual(s) to merge.
#' @param match The percent matching individuals must be before merging.
#' @param het The symbol(s) used for heterozygous calls.
#' 
#' @details
#' This function merges the calls from duplicate individuals. If a call for a given marker does not match between individuals, it is called heterozygous. All heterozygous markers are stored as "H".
#'
#' @return A gbs object with duplicate individuals merged.
#'
#' @examples
#' data(dh)
#' dh = hap.collapse(dh,names=c("TIGER","DANBY"),match=0.9)
#' head(dh$calls)
#' 
#' @export

hap.collapse <- function(hap, names, match, het="H") {
  
  if(class(hap)!="gbs") {
    stop("hap argument is incorrect type. Use hap.read to create gbs object.")
  }
  
  hap1 <- hap
  hap <- hap$calls
  
  if(missing(names)) {
    stop("Names of lines to merge must be specified.")
  } else {
    for(i in 1:length(unique(names))) {
      count <- length(colnames(hap)[grepl(names[i],colnames(hap),ignore.case = T)])
      
      if(count==0) {
        stop(paste("\"",names[i],"\""," not present in hap object.",sep=""))
      }
      
      if(count==1) {
        stop(paste("No duplicates for ","\"",names[i],"\""," found in hap object.",sep=""))
      }
    }
  }
  
  hap.orig <- hap
  
  # Removes lines that don't match each other at a certain percent identity
  if(!missing(match)) {
    
    hap.new <- matrix(nrow=nrow(hap),ncol=0)
    line.names <- names
    
    for(i in 1:length(unique(line.names))) {
      
      # Remove duplicates from original object
      hap.orig <- hap.orig[,!grepl(unique(line.names)[i],colnames(hap.orig),ignore.case = T)]
      
      hap.subset <- hap[,grepl(unique(line.names)[i],colnames(hap),ignore.case = T)]
      subset.identity <- allele.match(hap.subset,result = "percent")
      diag(subset.identity) <- NA
      subset.identity[lower.tri(subset.identity)] <- t(subset.identity)[lower.tri(subset.identity)]
      
      sub.array <- apply(subset.identity>=match,MARGIN = 1,function(x) length(x[x==TRUE]))
      sub.max <- max(apply(subset.identity>=match,MARGIN = 1,function(x) length(x[x==TRUE])))
      
      which(sub.array==sub.max)
      hap.subset <- hap.subset[,which(sub.array==sub.max)]
      
      if(length(which(sub.array==sub.max))==1) {
        print(paste("Only one individual from ",unique(line.names)[i], " matches at the specified threshold.",sep=""))
        hap.subset <- data.frame(hap.subset)
        names(hap.subset) <- unique(line.names)[i]
        hap.orig <- cbind(hap.subset,hap.orig)
        names <- names[names != unique(line.names)[i]]
      } else {
        hap.new <- cbind(hap.new,hap.subset)
      }
      
      print(paste(unique(line.names)[i],": ",length(which(sub.array==sub.max))," used out of ",length(sub.array),sep=""))
    }
    
    if(ncol(hap.new)==0) {
      stop("All individuals were removed.")
    }

  } else {
    match <- NULL
    hap.new <- hap
  }
  
  hap.out <- matrix(ncol=length(unique(names)),nrow=nrow(hap.new))
  
  for(i in 1:length(unique(names))) {
    
    # Remove all dupes from original hap object
    hap.orig <- hap.orig[,!grepl(unique(names)[i],colnames(hap.orig),ignore.case = T)]
    
    hap.subset <- t(hap.new[,grepl(unique(names)[i],colnames(hap.new),ignore.case = T)])
    x=apply(hap.subset,2,table)
    subset.calls <- lapply(x,call.fnc)
    hap.out[,i] <- unlist(subset.calls)
    print(paste(unique(names)[i],": ",nrow(hap.subset)," individuals were merged.", sep=""))
  }
  
  hap.out <- as.data.frame(hap.out)
  colnames(hap.out) <- unique(names)
  hap.orig <- cbind(hap.out,hap.orig)
  
  hap1$calls = hap.orig
  invisible(hap1)
}

call.fnc <- function(x) {
  count <- sum(x)
  if(length(x) == 1) {
    return(names(x[1]))
  }
  if(length(x) == 2 && "N"%in%names(x)) {
    return(names(x)[names(x)!="N"])
  }
  if("H"%in%names(x)) { # TODO change for IUPAC
    return("H")
  }
  return("H")
}