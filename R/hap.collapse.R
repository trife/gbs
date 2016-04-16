#' Hap Collapse
#' 
#' Combines multiple genotypes with the same name in a hap object
#' 
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' 
#' @param hap the input hap object
#' @param names the name of a specific line to merge
#' @param match percent identical genotypes must be before merging
#' @param method method used to determine which lines to remove when they don't match each other
#' 
#' @keywords
#' 
#' @examples
#' 
#' @export

hap.collapse <- function(hap,names,match,method){
  
  #TODO make call.fnc more robust
  #TODO add option for matching percent, call allele.match to calculate
  #TODO attempts to automatically identify names to merge will break if names contain periods
  #TODO add different methods for removing lines (strict, step, etc.)
  
  call.fnc = function(x) {
    count = sum(x)
    if(length(x) == 1) {
      return(names(x[1]))
    }
    if(length(x) == 2 && "N"%in%names(x)) {
      return(names(x)[names(x)!="N"])
    }
    if("H"%in%names(x)) {
      return("H")
    }
    return("H")
  }
  
  if(missing(names)) {
    # Substring names by period in case columns were originally identically named
    names = unlist(lapply(strsplit(colnames(hap),".",fixed=TRUE), `[[`, 1))
    
    if(!any(duplicated(names))) {
      stop("No duplicate names found in the hap object.")
    }
  }
  
  # Removes lines that don't match each other at a certain percent identity
  if(!missing(match)) {
    
    hap.new = matrix(nrow=nrow(hap),ncol=0)
    
    for(i in 1:length(unique(names))) {
      hap.subset = hap[,grepl(unique(names)[i],colnames(hap))]
      subset.identity = allele.match(hap.subset,result = "percent")
      diag(subset.identity) <- NA
      subset.identity[lower.tri(subset.identity)] <- t(subset.identity)[lower.tri(subset.identity)]
      
      sub.array = apply(subset.identity>=match,MARGIN = 1,function(x) length(x[x==TRUE]))
      sub.max = max(apply(subset.identity>=match,MARGIN = 1,function(x) length(x[x==TRUE])))
      
      which(sub.array==sub.max)
      hap.subset = hap.subset[,which(sub.array==sub.max)]
      hap.new = cbind(hap.new,hap.subset)
      
      print(paste(unique(names)[i],": ",length(which(sub.array==sub.max))," used out of ",length(sub.array),sep=""))
    }
    
    if(ncol(hap.new)==0) {
      stop("All individuals were removed.")
    }

  } else {
    match <- NULL
    hap.new = hap
  }
  
  hap.out = matrix(ncol=length(unique(names)),nrow=nrow(hap.new))
  
  for(i in 1:length(unique(names))) {
    hap.subset = t(hap.new[,grepl(unique(names)[i],colnames(hap.new))])
    x=apply(hap.subset,2,table)
    subset.calls = lapply(x,call.fnc)
    hap.out[,i] = unlist(subset.calls)
    print(paste(unique(names)[i],": ",nrow(hap.subset)," individuals were merged.",sep=""))
  }
  
  # TODO this currently deletes everything else; needs to remove columns from the original dataset and return it instead
  
  hap.out = as.data.frame(hap.out)
  colnames(hap.out) = unique(names)

  hap.out
}