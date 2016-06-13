#' summarize the gbs object.
#' 
#' Processes the gbs object and calculates allele statistics. Returns allele stats, a geno object, and graphs statistics.
#' 
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' @author Jesse Poland, \email{jpoland@@ksu.edu}
#' 
#' @param hap The gbs object to be processed.
#' @param geno A logical value that will convert the marker calls to a numeric geno format.
#' @param encoding The numbers used to encode alleles (minor, het, major).
#' @param graph A logical value that will output graphs for summary statistics (blank wells, maf, percent present, het, percent het).
#' @param het The symbol(s) used for heterozygous calls.
#' @param missing The symbol(s) used for missing calls.
#' 
#' @details
#' This function checks for and removes any wells that contain "blank" in the sample name, calculates allele counts, percent present, and minor allele frequency, and optionally creates a numeric geno object (stored as hap$geno). The encoding argument can be used to change how the minor, major, and heterozygous alleles are encoded. The het and missing arguments can be a character or vector to allow for multiple matches (e.g. IUPAC het = c("R","Y","S","W","K","M")).
#'
#' @examples
#' data(wheat)
#' hap = gbs.summary(hap,geno=T,graph=T)
#'
#' @export

gbs.summary <- function(hap, geno=F, encoding=c(-1,0,1), graph=F, het="H", missing="N") {
  
  if(class(hap)!="gbs") {
    stop("hap argument is incorrect type. Use hap.read to create gbs object.")
  }
  
  # Check blank wells for data
  if(any(grepl("blank",colnames(hap$calls),ignore.case=TRUE))) {
    missing.blank <- hap$calls[,grepl("BLANK",colnames(hap$calls),ignore.case=TRUE)]=="N"
    blank <- as.matrix(apply(!missing.blank, 2, sum))
    print("Number of reads in blank wells: ")
    print(blank)
    
    if(graph) {
      missing.blank <- hap$calls[,grepl("BLANK",colnames(hap$calls),ignore.case=TRUE)]=="N"
      blank <- as.matrix(apply(!missing.blank, 2, sum))
      snptot <- colSums(hap$calls!="N")
      upper_limit <- round(max(snptot)/1000)*1000
      hist(snptot, xlim=c(0,upper_limit), breaks=seq(0,upper_limit, by=1000), main="SNPs number per Sample", xlab="Number of SNPs", sub="Blank wells in red")
      hist(blank, col="red",xlim=c(0,upper_limit), breaks=seq(0,upper_limit, by=1000),add=TRUE )
      hap$calls <- hap$calls[,!grepl("blank",colnames(hap$calls), ignore.case=TRUE)]
    }
    
    print("Removing blank wells...")
    
    # Remove blank columns
    hap$calls <- hap$calls[,!grepl("blank",colnames(hap$calls), ignore.case=TRUE)]
  }
  
  # Calculate allele stats
  print("Calculating statistics...")
  
  a <- substring(hap$header$alleles,1,1)
  b <- substring(hap$header$alleles,3,3)
  
  counts <- apply(hap$calls,MARGIN=1,function(x) table(t(x)))
  
  alleleA <- lapply(1:length(counts),function(x) counts[[x]][names(counts[[x]])==a[x]])
  alleleB <- lapply(1:length(counts),function(x) counts[[x]][names(counts[[x]])==b[x]])
  
  calc.het <- lapply(1:length(counts),function(x) counts[[x]][names(counts[[x]]) %in% het])
  calc.missing <- lapply(1:length(counts),function(x) counts[[x]][names(counts[[x]]) %in% missing])
  
  alleleA[!(sapply(alleleA, length))] <- 0
  alleleB[!(sapply(alleleB, length))] <- 0
  calc.het[!(sapply(calc.het, length))] <- 0
  calc.missing[!(sapply(calc.missing, length))] <- 0
  
  alleleA <- unlist(alleleA)
  alleleB <- unlist(alleleB)
  calc.het <- unlist(calc.het)
  calc.missing <- unlist(calc.missing)
  
  present <- (alleleA + alleleB) / (alleleA+alleleB+calc.missing+calc.het)
  
  # Calculate and add summary stats
  maf <- apply(cbind(alleleA, alleleB), 1, min) / apply(cbind(alleleA, alleleB, calc.het), 1, sum)
  phet <- calc.het / apply(cbind(alleleA, alleleB, calc.het), 1, sum)
  present <- as.numeric(as.character(present))
  hap$stats <- data.frame(alleleA,alleleB,calc.het,calc.missing,present,maf,phet)
  
  if(geno==TRUE) {
    print("Converting to geno...")
    
    ## Define allele a and allele b
    a[hap$stats$alleleA<hap$stats$alleleB] <- substring(hap$header$alleles,3,3)[hap$stats$alleleA<hap$stats$alleleB]
    b[hap$stats$alleleA<hap$stats$alleleB] <- substring(hap$header$alleles,1,1)[hap$stats$alleleA<hap$stats$alleleB]
    
    # Default: turn minor into -1, major into 1, het into 0
    hap01 <- hap$calls
    hap01[,1:ncol(hap01)] <- NA
    
    print("Converting major allele...")
    hap01[hap$calls == a] <- encoding[1]
    
    print("Converting minor allele...")
    hap01[hap$calls == b] <- encoding[3]
    
    print("Converting het allele...")
    hap01[hap$calls %in% het] <- encoding[2]
    
    #Convert to geno
    geno <- as.matrix(hap01)
    class(geno) <- "numeric"
    geno <- t(geno)
    hap$geno <- geno
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
    
    if("calc.het"%in%colnames(hap$stats)) {
      hist(calc.het, main="Number of Heterozygotes", xlab="Number of heterozygous per SNP loci", ylab="Number of SNPs")
    }
    
    if("phet"%in%colnames(hap$stats)) {
      hist(phet, main="Percent Heterozygous", xlab="Percent Heterozygous", ylab="Number of SNPs")
    }
    
    par(mfrow=c(1,1))
  }
  
  invisible(hap)
}