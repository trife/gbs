#' Filter and summarize the hap object.
#' 
#' Processes the hap object, calculates summary statistics, returns list with hap and/or geno objects, and graphs several summary statistics.
#' 
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' @author Jesse Poland, \email{jpoland@@ksu.edu}
#' 
#' @param hap The gbs object to be processed.
#' @param geno A logical value that will convert the marker calls to a numeric geno format.
#' @param graph A logical value that will output graphs for summary statistics (blank wells, maf, percent present, het, percent het).
#' 
#' @details
#'
#' @examples
#' data(wheat)
#' hap = hap.read(wheat)
#' hap = filter.summary(hap,project="wheat",output="hap",graph=T)
#'
#' @export

filter.summary <- function(hap, geno=F, graph=F){
  
  # TODO integrate IUPAC https://bytebucket.org/tasseladmin/tassel-5-source/wiki/docs/Tassel5UserGuide.pdf
  # TODO use proper hap format https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load
  
  if(class(hap)!="gbs") {
    stop("hap object is incorrect type. Use hap.read to create gbs object.")
  }
  
  # Check blank wells for data
  if(any(grepl("blank",colnames(hap$calls),ignore.case=TRUE))) {
    missing.blank = hap$calls[,grepl("BLANK",colnames(hap$calls),ignore.case=TRUE)]=="N"
    blank = as.matrix(apply(!missing.blank, 2, sum))
    print("Number of reads in blank wells: ")
    print(blank)
    
    if(graph) {
      missing.blank = hap$calls[,grepl("BLANK",colnames(hap$calls),ignore.case=TRUE)]=="N"
      blank = as.matrix(apply(!missing.blank, 2, sum))
      snptot=colSums(hap$calls[,data.col:ncol(hap$calls)]!="N")
      upper_limit = round(max(snptot)/1000)*1000
      hist(snptot, xlim=c(0,upper_limit), breaks=seq(0,upper_limit, by=1000), main="SNPs number per Sample", xlab="Number of SNPs", sub="Blank wells in red")
      hist(blank, col="red",xlim=c(0,upper_limit), breaks=seq(0,upper_limit, by=1000),add=TRUE )
      hap$calls = hap$calls[,!grepl("blank",colnames(hap$calls), ignore.case=TRUE)]
    }
    
    print("Removing blank wells...")
    
    # Remove blank columns
    hap$calls = hap$calls[,!grepl("blank",colnames(hap$calls), ignore.case=TRUE)]
  }
  
  # Calculate allele stats
  a = substring(hap$header$alleles,1,1)
  b = substring(hap$header$alleles,3,3)
  
  counts = apply(hap$calls,MARGIN=1,function(x) table(t(x)))
  
  alleleA = lapply(1:length(counts),function(x) counts[[x]][names(counts[[x]])==a[x]])
  alleleB = lapply(1:length(counts),function(x) counts[[x]][names(counts[[x]])==b[x]])
  
  het = lapply(1:length(counts),function(x) counts[[x]][names(counts[[x]])=="H"]) #TODO change for IUPAC
  missing = lapply(1:length(counts),function(x) counts[[x]][names(counts[[x]])=="N"]) #TODO change for IUPAC
  
  alleleA[!(sapply(alleleA, length))] = 0
  alleleB[!(sapply(alleleB, length))] = 0
  het[!(sapply(het, length))] = 0
  missing[!(sapply(missing, length))] = 0
  
  alleleA = unlist(alleleA)
  alleleB = unlist(alleleB)
  het = unlist(het)
  missing = unlist(missing)
  
  present = (alleleA + alleleB) / (alleleA+alleleB+missing+het)
  
  # Calculate and add summary stats
  maf = apply(cbind(alleleA, alleleB), 1, min) / apply(cbind(alleleA, alleleB, het), 1, sum)
  phet = het / apply(cbind(alleleA, alleleB, het), 1, sum)
  present = as.numeric(as.character(present))
  hap$stats = data.frame(alleleA,alleleB,het,missing,present,maf,phet)
  
  if(geno==TRUE) {
    
    ## Define allele a and allele b
    a[hap$stats$alleleA<hap$stats$alleleB] = substring(hap$header$alleles,3,3)[hap$stats$alleleA<hap$stats$alleleB]
    b[hap$stats$alleleA<hap$stats$alleleB] = substring(hap$header$alleles,1,1)[hap$stats$alleleA<hap$stats$alleleB]
    
    ## Turn allele a and allele b into -1 and 1.  Het into 0
    hap01 = hap$calls
    hap01 = NA
    hap01[hap$calls == a] = "-1"
    hap01[hap$calls == b] = "1"
    hap01[hap$calls == "H"] = "0" # TODO change for IUPAC
    
    #Convert to geno
    geno = as.matrix(hap01)
    class(geno) = "numeric"
    geno = t(geno)
    hap$geno = geno
  }
  
  if(graph) {
    par(mfrow=c(2,2))
    
    # Graph population parameters
    if("maf"%in%colnames(hap$stats)) {
      hist(maf, main="Minor Allele Frequency", xlab="MAF Value", ylab="Number of SNPs")
    }
    
    if("present"%in%colnames(hap$stats)) {
      hist(present, main="% Present of Each SNP", xlab="Percent Present", ylab="Number of SNPs")
    }
    
    if("het"%in%colnames(hap$stats)) {
      hist(het, main="Number of Heterozygotes", xlab="Number of heterozygous per SNP loci", ylab="Number of SNPs")
    }
    
    if("phet"%in%colnames(hap$stats)) {
      hist(phet, main="Percent Heterozygous", xlab="Percent Heterozygous", ylab="Number of SNPs")
    }
    
    par(mfrow=c(1,1))
  }
  
  invisible(hap)
}