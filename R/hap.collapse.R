#' Hap Collapse
#' 
#' Combines multiple genotypes with the same name in a hap object
#' 
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' 
#' @param hap the input hap object
#' 
#' @keywords
#' 
#' @examples
#' 
#' @export

hap.collapse <- function(hap){

  if(!any(duplicated(colnames(hap)))) {
    stop("No common individual names found in the hap object.")
  }

}