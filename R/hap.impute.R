#' Impute markers in gbs
#' 
#' Takes in hap object and impute markers 
#' 
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' 
#' @param hap.obj the hap object to impute
#' @param method the method to use for imputation
#' @param parents columns for the two parents being used to convert to RQTL and AB format
#' 
#' @keywords 
#' 
#' @examples
#' 
#' @export

hap.impute <- function(hap,method=c("mean","EM","RF"),parents=null){
  require(rrBLUP)
  hap.obj = hap

  mean.F = function(...){
    
  }
  
  EM.F = function(...){
    
  }
  
  RF.F = function(...){
    
  }
  
  if(any(methods=="mean")){
    output$MEGA4 = RQTL.F()
  }
  
  if(any(methods=="EM")){
    output$MEGA4 = AB.F()
  }
  
  if(any(methods=="RF")){
    output$MEGA4 = GAPIT.F()
  }
  
  #TODO change output types
}