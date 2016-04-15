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
    if(length(parents) != 2) {
      stop("Exactly two parents must be specified for AB format.")
    }
    
    # Check for multiple occurrences of parents, and that they both exist
    if(is.character(parents)) {
      if(length(grep(paste(parents,collapse="|"), colnames(hap), value=TRUE))<2) {
        stop("Unable to find both parent data columns.")
      }
      
      if(length(grep(paste(parents,collapse="|"), colnames(hap), value=TRUE))>2) {
        stop("Found too many parent data columns. Use hap.collapse to merge.")
      }
    }

    # If user specifies column number instead of names
    if(is.numeric(parents)) {
      p1col = parents[1]
      p2col = parents[2]
      
      print(paste("Using ",colnames(hap)[parents[1]]," and ",colnames(hap)[parents[2]]," as parents.",sep=""))
    } else {
      p1col = which(grepl(parents[1],colnames(hap)))
      p2col = which(grepl(parents[2],colnames(hap)))
    }
    
    ## Remove markers that are hets, identical, or missing in both parents
    print("Removing markers that are het, identical, or missing in both parents...")

    hap.clean = hap[hap[,p1col]!="H" & hap[,p2col]!="H",]
    hap.clean = hap.clean[hap.clean[,p1col]!="N" & hap.clean[,p2col]!="N",]
    hap.clean = hap.clean[(hap.clean[,p1col]!=hap.clean[,p2col]),]

    print(paste("Removed ",nrow(hap)-nrow(hap.clean)," markers.",sep=""))
    
    ## Identify two alleles for each marker
    alleles = apply(hap.clean,MARGIN = 1,function(x) names(table(unlist(lapply(x,as.character)))))
    alleles = lapply(alleles,function(x) x[!x%in%c("H","N")]) # TODO add option for ignorable genotypes
    
    ## Identify which allele belongs to which parent and impute missing from the other parent
    # TODO add an if here to skip this if there's no missing data
    p1 = hap.clean[,p1col]
    p2 = hap.clean[,p2col]
    
    allele1 = lapply(alleles, `[[`, 1)
    allele1 = unlist(allele1)
    allele2 = lapply(alleles, `[[`, 2)
    allele2 = unlist(allele2)
    
    hap.ab = hap.clean
    
    hap.ab[,p1col] = as.character(hap.ab[,p1col])
    hap.ab[,p2col] = as.character(hap.ab[,p2col])
    
    print("Imputing parental genotypes when only one present...")
    
    for(i in 1:length(alleles)){
      if (hap.ab[,p1col][i] == "N") {
        if (hap.ab[,p2col][i] == allele1[i]) {
          hap.ab[,p1col][i] = allele2[i]
        } else if (hap.ab[,p2col][i] == allele2[i]) {
          hap.ab[,p2col][i] = allele1[i]
        }
      } else if (hap.ab[,p2col][i] == "N") {
        if (hap.ab[,p1col][i] == allele1[i]) {
          hap.ab[,p2col][i] = allele2[i]
        } else if (hap.ab[,p1col][i] == allele2[i]) {
          hap.ab[,p2col][i] = allele1[i]
        }
      }
    }
    
    ## Convert to A/B (converted to Y/Z since A is also a nucleotide)
    hap.calls = hap.ab
    hap.calls <- data.frame(lapply(hap.calls, as.character), stringsAsFactors=FALSE, check.names=FALSE)
    hap.calls[hap.calls==hap.ab[,p1col]] = "Y"
    hap.calls[hap.calls==hap.ab[,p2col]] = "Z"
    hap.calls[hap.calls=="Y"] = "A"
    hap.calls[hap.calls=="Z"] = "B"
    
    hap.calls
  }
  
  GAPIT.F = function(...) {
    
  }
  
  output = list()
  
  if(any(format=="MEGA4")){
    output[[MEGA4]] = MEGA4.F()
  }
  
  if(any(format=="FST")){
    output$FST = FST.F()
  }
  
  if(any(format=="RQTL")){
    output$RQTL = RQTL.F()
  }
  
  if(any(format=="AB")){
    output$AB = AB.F()
  }
  
  if(any(format=="GAPIT")){
    output$GAPIT = GAPIT.F()
  }
  
  names(output) = format
  cat("DONE","\n")
  invisible(output)
}