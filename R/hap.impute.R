#' Impute markers in gbs
#' 
#' Takes in hap object and impute markers 
#' 
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' 
#' @param hap.obj the hap object to convert
#' @param delim the deliminator used in the file
#' @param exp.format the format you wish to convert the hap object to
#' @param parents columns for the two parents being used to convert to RQTL and AB format
#' 
#' @keywords 
#' 
#' @examples
#' 
#' @export

hap.impute <- function(hap,delim,exp.format=c("MEGA4","FST","STRUCTURE","RQTL","AB","GAPIT"),parents=null){
  hap.obj = hap

  MEGA4.F = function(...){
    
  }
  
  FST.F = function(...){
    
  }
  
  RQTL.F = function(...){
    
  }
  
  AB.F = function(...){
    
  }
  
  GAPIT.F = function(...){
    
  }
  
  if(any(methods=="MEGA4")){
    output$MEGA4 = MEGA4.F()
  }
  
  if(any(methods=="FST")){
    output$MEGA4 = FST.F()
  }
  
  if(any(methods=="RQTL")){
    output$MEGA4 = RQTL.F()
  }
  
  if(any(methods=="AB")){
    output$MEGA4 = AB.F()
  }
  
  if(any(methods=="GAPIT")){
    output$MEGA4 = GAPIT.F()
  }
  
  #TODO change output types
}