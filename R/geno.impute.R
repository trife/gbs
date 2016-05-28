#' Impute GBS markers.
#' 
#' Takes in a geno object and imputes markers using various methods.
#' 
#' @author Chris Gaynor, \email{chris.gaynor@roslin.ed.ac.uk}
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' @author Jessica Rutkoski, \email{jer263@cornell.edu}
#' 
#' @param geno The geno object to impute. Each row is an individual and each column is a marker. Genotypic data is expected to be encoded as \code{-1,0,1}.
#' @param method The method to use for imputation. See details.
#' @param k The number of clusters for KNN imputation.
#' @param maxiter The maximum number of iterations for RF and SVD imputation.
#' @param n.core The number of cores to use for RF imputation.
#' @param ... Additional arguments for rrBLUP::A.mat, randomForest::na.roughfix, missForest::missForest
#' 
#' @details
#' @section Methods:
#' Markers can be imputed using mean, median, EM, Random Forest (RF), K Nearest Neighbors (KNN), and Singular Value Decomposition (SVD)
#' 
#' @return
#'
#' @examples
#' data(wheat)
#' hap = hap.read(wheat)
#' geno = filter.summary(hap,project="wheat",output=c("geno","hap"))$geno
#' geno.impute = hap.impute(geno,c("MEAN","EM"))
#' head(geno.impute$mean)
#' head(geno.impute$EM)
#'
#' @export

geno.impute <- function(geno, method, k=4, maxiter=10, n.core=1, ...){
  
  method = toupper(method)
  
  all.methods = c("MEAN","EM","RF","MEDIAN","KNN","SVD")
  
  if(!any(method%in%all.methods) || missing(method)) {
    stop("Method missing or incorrectly specified.")
  }
  
  if(!any(apply(geno,MARGIN=2,is.numeric()))) {
    stop("Geno object contains non-numeric elements.")
  }
  
  if(.Platform$OS.type=="windows") { 
    n.core=1
    parallelize = "no"
  } else {
    parallelize = "variables"
  }
  
  data = list()

  if(any(method=="MEAN")) {
    print("Imputing markers using mean...")
    data$mean = rrBLUP::A.mat(geno, impute.method="mean", n.core=n.core, return.imputed=TRUE, ...)$imputed
  }
  
  if(any(method=="MEDIAN")) {
    print("Imputing markers using median...")
    data$median = randomForest::na.roughfix(geno)
  }
  
  if(any(method=="EM")) {
    print("Imputing markers using EM...")
    data$EM = rrBLUP::A.mat(geno, impute.method="EM", n.core=n.core, return.imputed=TRUE, ...)$imputed
  }
  
  if(any(method=="RF")) {
    print("Imputing markers using RF...")
    cl = parallel::makeCluster(n.core)
    doParallel::registerDoParallel(cl)
    data$RF = missForest::missForest(geno, parallelize=parallelize, maxiter=maxiter, ...)$ximp
    parallel::stopCluster(cl)
  }
  
  if(any(method=="KNN")) {
    print("Imputing markers using KNN...")
    data$KNN = kNNI(x=geno, k=k)$x
  }

  if(any(method=="SVD")) {
    print("Imputing markers using SVD...")
    data$SVD = bcv::impute.svd(x=geno, maxiter=maxiter, ...)$x
  }
  
  invisible(data)
}

kNNI <- function (x, k, verbose = F, x.dist=NULL) {
  x= t(x)
  indexkeep=1:nrow(x) 
  rg= range(na.omit(as.vector(x)))
  recode= function(vec, rg){
    ind2=grep(TRUE, vec==rg[2])
    ind1=grep(FALSE, vec==rg[2])
    vec[ind2]=rg[1]
    vec[ind1]=rg[2]
    return(vec)
  }
  x= rbind(x, apply(x, 2, recode, rg=rg) )
  
  
  if (k >= nrow(x)) 
    stop("k must be less than the number of rows in x")
  missing.matrix = is.na(x)
  numMissing = sum(missing.matrix)
  if (verbose) {
    print(paste("imputing on", numMissing, "missing values with matrix size", 
                nrow(x) * ncol(x), sep = " "))
  }
  if (numMissing == 0) {
    return(x)
  }
  
  if ( is.matrix(x.dist)==FALSE ){ 
    if (verbose) 
      print("Computing distance matrix...")
    x.dist = as.matrix(dist(x, upper = T))
    if (verbose) 
      print("Distance matrix complete")
  }
  
  missing.rows.indices = which(apply(missing.matrix, 1, function(i) {
    any(i)
  }))
  x.missing = (cbind(1:nrow(x), x))[missing.rows.indices, ]
  x.missing.imputed = t(apply(x.missing, 1, function(i) {
    rowIndex = i[1]
    i.original = i[-1]
    if (verbose) 
      print(paste("Imputing marker", rowIndex, sep = " "))
    missing.cols = which(missing.matrix[rowIndex, ])
    if (length(missing.cols) == ncol(x)) 
      warning(paste("Marker", rowIndex, "is completely missing", 
                    sep = " "))
    imputed.values = sapply(missing.cols, function(j) {
      neighbor.indices = which(!missing.matrix[, j])
      knn.ranks = order(x.dist[rowIndex, neighbor.indices])
      
      if (x.dist[rowIndex, neighbor.indices][knn.ranks][1:k][1]==0){
        knn = neighbor.indices[(knn.ranks[1])]
        x[knn, j]
      }else{
        
        wghts= 1/(x.dist[rowIndex, neighbor.indices][knn.ranks][1:k])^2
        wghts= wghts/sum(wghts)
        knn = neighbor.indices[(knn.ranks[1:k])]
        sum(x[knn, j]*wghts )
      }
    })
    i.original[missing.cols] = imputed.values
    i.original
  }))
  x[missing.rows.indices, ] = x.missing.imputed
  missing.matrix2 = is.na(x)
  x[missing.matrix2] = 0
  x= t(x)
  x= x[,indexkeep]
  return(list(x = x, missing.matrix = missing.matrix))
}