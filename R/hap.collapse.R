#' Hap Collapse
#' 
#' Combines multiple genotypes with the same name in a hap object
#' 
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' 
#' @param hap the input hap object
#' @param match percent identical genotypes must be before merging
#' @param names the name of a specific line to merge
#' 
#' @keywords
#' 
#' @examples
#' 
#' @export

hap.collapse <- function(hap,match,names){
  
  #TODO make call.fnc more robust
  #TODO add option for matching percent, call allele.match to calculate
  #TODO attempts to automatically identify names to merge will break if names contain periouds
  
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
  
  # TODO Remove lines that don't match each other at a certain percent identity
  if(!missing(match)) {
    for(i in 1:length(unique(names))) {
      hap.subset = hap[,grepl(unique(names)[i],colnames(hap))]
      subset.identity = allele.match(hap.subset,result = "percent")
      diag(subset.identity) = NA
    }
  } else {
    match <- NULL
  }
  
  hap.out = matrix(ncol=length(unique(names)),nrow=nrow(hap))
  
  for(i in 1:length(unique(names))) {
    hap.subset = t(hap[,grepl(unique(names)[i],colnames(hap))])
    x=apply(hap.subset,2,table)
    subset.calls = lapply(x,call.fnc)
    hap.out[,i] = unlist(subset.calls)
    print(paste(unique(names)[i],": ",nrow(hap.subset)," individuals were merged.",sep=""))
  }
  
  hap.out = as.data.frame(hap.out)
  colnames(hap.out) = unique(names)

  hap.out
}