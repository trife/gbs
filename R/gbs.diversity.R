#' GBS Diversity
#' 
#' Calculate different measures of diveristy in a hap object
#' 
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' 
#' @param hap the hap object to manipulate
#' @param delim the delimiter for the hap object (if file)
#' 
#' @keywords
#' 
#' @examples
#' 
#' @export

gbs.diversity <- function(geno){
  require(rrBLUP)
  
  
  ## TODO diversity, etc.
  ae = 1/(1-hap$percentHET)
  ae2 = 1/(hap$MAF^2)
  nei = (1-hap$MAF^2-(1-hap$MAF)^2)
}