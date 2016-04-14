#' Impute markers in gbs
#' 
#' Takes in hap object and impute markers 
#' 
#' @author Chris Gaynor
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

hap.impute <- function(hap,method=c("mean","EM","RF","median","KNN","SVD"),parents=null){
  hap.obj = hap

  if(any(methods=="mean")){
    data[[1]] = rrBLUP::A.mat(hap.obj, impute.method="mean", n.core=n.core, return.imputed=TRUE, ...)$imputed
  }
  
  if(any(methods=="median")) {
    data[[1]] = randomForest::na.roughfix(hap.obj)
  }
  
  if(any(methods=="EM")){
    data[[1]] = rrBLUP::A.mat(hap.obj, impute.method="EM", n.core=n.core, return.imputed=TRUE, ...)$imputed
  }
  
  if(any(methods=="RF")){
    cl = parallel::makeCluster(n.core)
    registerDoParallel(cl)
    data[[1]] = missForest::missForest(hap.obj, parallelize="variables", ...)$ximp
    parallel::stopCluster(cl)
    rm(cl)
  }
  
  if(any(methods=="KNN")){
    #TODO
  }
  
  if(any(methods=="SVD")){
    #TODO
  }
  
  invisible(data)
}