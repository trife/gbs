#' Convert hap objects to other formats
#' 
#' Takes in hap files and converts them to the formats needed for other analysis programs 
#' 
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#' 
#' @param hap.obj the hap object to convert
#' @param format the format you wish to convert the hap object to
#' @param parents columns for the two parents being used to convert to RQTL and AB format
#' 
#' @keywords 
#' 
#' @examples
#' 
#' @export

hap.convert <- function(hap,format=c("MEGA4","FST","STRUCTURE","RQTL","AB","GAPIT"),parents=NULL) {
  hap.obj = hap
  
  if(missing(format)) {
    stop("No export format specified.")
  }

  MEGA4.F = function(...) {
    
  }
  
  FST.F = function(...) {
    
  }
  
  RQTL.F = function(...) {
    
  }
  
  AB.F = function(...) {
    if(missing(parents) | length(parents) = 1) {
      stop("Two parents must be specified for AB format.")
    }
    
    # Try to use parents as column numbers in hap matrix else find names
    if(is.numeric(parents)) {
      
    } else {
      
    }
  }
  
  GAPIT.F = function(...) {
    
  }
  
  if(any(exp.format=="MEGA4")){
    output$MEGA4 = MEGA4.F()
  }
  
  if(any(exp.format=="FST")){
    output$MEGA4 = FST.F()
  }
  
  if(any(exp.format=="RQTL")){
    output$MEGA4 = RQTL.F()
  }
  
  if(any(exp.format=="AB")){
    output$MEGA4 = AB.F()
  }
  
  if(any(exp.format=="GAPIT")){
    output$MEGA4 = GAPIT.F()
  }
  
  names(output) = format
  cat("DONE","\n")
  invisible(output)
}