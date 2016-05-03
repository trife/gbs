#' Read hap file(s) or hap object.
#'
#' Read a file or data frame to create a hap object
#'
#' @author Trevor Rife, \email{trife@@ksu.edu}
#'
#' @param hap The hap object, file, or folder containing hap files.
#' @param delim The delimter used in the file(s).
#' @param data.col The column number for the first individual. This removes columns not required if specified.
#' @param col.names The required columns for the gbs package to function correctly. See \code{Details}.
#' @param genotypes A vector of marker calls that are used in the hap object.
#'
#' @details
#' If a folder is passed to \code{hap.read}, \code{hap.join} is called to merge all files within the folder.
#' 
#' \code{col.names} specifies the names of the columns corresponding to the tag, the alleles (e.g. A/G format), the position of the marker in the tag, the physical or genetic position, and the chromosome, in that order. All of these columns are required, even if blank. If columns are named differently, pass the names via the col.names option.
#'
#' \code{genotypes} allows for non-standard genotype calls to be included in the file. A basic data integrity checks to ensure that the resulting data frame only contains genotypic calls.
#'
#' @return A verified hap data frame.
#'
#' @examples
#' data(wheat)
#' hap = hap.read(wheat)
#'
#' @export

hap.read <- function(hap.obj, delim="\t", data.col, col.names = c("rs","alleles","tagpos","pos","chrom"), genotypes=c(NA,"A","T","C","G","H","N")){
  
  if(class(hap.obj)!= "data.frame") {
    if(!file.exists(hap.obj)) {
      stop("File or folder does not exist.")
    }
    
    if(file.info(hap.obj)$isdir) {
      hap = hap.join(hap.obj,delim)
    }
    
    if(!file.info(hap.obj)$isdir) {
      hap = data.table::fread(input=hap.obj, sep=delim,check.names=FALSE,header=TRUE, data.table = F,strip.white=T)
      
      if(!"dif"%in%colnames(hap)) {
        hap = cbind(dif=1,hap)  
      }
    }
  } else {
    hap = hap.obj
    
    if(!"dif"%in%colnames(hap)) {
      hap = cbind(dif=1,hap)  
    }
  }

  if(!col.names[1]%in%colnames(hap)) {
    stop("Tag column (rs) missing from hap object")
  }
  
  if(!col.names[2]%in%colnames(hap)) {
    stop("Alleles column (alleles) missing from hap object")
  }
  
  if(!col.names[3]%in%colnames(hap)) {
    stop("Marker position column (tagpos) missing from hap object")
  }
  
  if(!col.names[4]%in%colnames(hap)) {
    stop("Physical position column (pos) missing from hap object")
  }
  
  if(!col.names[5]%in%colnames(hap)) {
    stop("Chromosome column (chrom) missing from hap object")
  }
  
  if(!missing(data.col)) {
    hap = cbind(dif=hap$dif,rs=hap$rs,alleles=hap$alleles,pos=hap$pos,chrom=hap$chrom,hap[,(data.col+1):ncol(hap)])
    
    # Simple data check for non-genotypes in first 10 roww
    if(!all(apply(hap[1:10,6:ncol(hap)],MARGIN=2,function(x) x%in%genotypes))) {
      stop("Non genotypes detected in hap object. Edit the genotypes parameter or check your hap object.")
    }
  }
  
  # Simple data check for non-genotypes in first 10 rows
  if(!all(apply(hap[1:10,6:ncol(hap)],MARGIN=2,function(x) x%in%genotypes))) {
    stop("Non genotypes detected in hap object. Edit the genotypes parameter or check your hap object.")
  }
  
  invisible(hap)
}