#' GBS Graphs
#' 
#' Calculates and graphs relevant stats for a hap object
#' 
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' @author Jared Crain, \email{jcrain@@ksu.edu}
#' 
#' @param hap the hap object to be processed
#' @param geno the geno object to be processed
#' @param file the output file name
#' @param data.col the column where data starts
#' 
#' @keywords
#' 
#' @examples
#' 
#' @export


gbs.graph <- function(hap,geno,file,data.col=12) {
  
  ## TODO scale page sizes by the number of lines or something
  if(!missing(file)) {
    pdf(file)
  }
  
  # Check for blank wells, make histogram if they exist
  if(any(grepl("BLANK",colnames(hap),ignore.case=TRUE))) {
    missing.blank = hap[,grepl("BLANK",colnames(hap),ignore.case=TRUE)]=="N"
    blank = as.matrix(apply(!missing.blank, 2, sum))
    snptot=colSums(hap[,data.col:ncol(hap)]!="N")
    upper_limit = round(max(snptot)/1000)*1000
    
    hist(snptot, xlim=c(0,upper_limit), breaks=seq(0,upper_limit, by=1000), main="SNPs number per Sample", xlab="Number of SNPs", sub="Blank wells in red")
    hist(blank, col="red",xlim=c(0,upper_limit), breaks=seq(0,upper_limit, by=1000),add=TRUE )
    
    hap = hap[,!grepl("blank",colnames(hap), ignore.case=TRUE)]
  }
  
  # Graph population parameters
  if(!"maf"%in%colnames(hap)) {
    print("No minor allele frequency (maf) column in hap object.")
  } else {
    hist(hap$maf, main="Minor Allele Frequency", xlab="MAF Value", ylab="Number of SNPs")
  }
  
  if(!"present"%in%colnames(hap)) {
    print("No present column in hap object.")
  } else {
    hist(hap$present, main="% Present of Each SNP", xlab="Percent Present", ylab="Number of SNPs")
  }
  
  if(!"het"%in%colnames(hap)) {
    print("No het column in hap object.")
  } else {
    hist(hap$het, main="Number of Heterozygotes", xlab="Number of heterozygous per SNP loci", ylab="Number of SNPs")
  }
  
  if(!"phet"%in%colnames(hap)) {
    print("No percent het (phet) column in hap object.")
  } else {
    hist(hap$phet, main="Percent Heterozygous", xlab="Percent Heterozygous", ylab="Number of SNPs")
  }
  
  # Make simple dendrogram
  if(!missing(geno)) {
    hc <- stats::hclust(dist(geno))
    plot(hc, main="GBS Dendroram",cex=.5)
  }
  
  # Make level plot of all filtered data
  fmatch=allele.match(hap[,data.col:ncol(hap)],result = "percent")
  fmatch[lower.tri(fmatch)] <- t(fmatch)[lower.tri(fmatch)]
  stats::heatmap(fmatch,cexRow=.5,cexCol=.5, main = "GBS Heatmap")

  if(!missing(file)) {
    dev.off()
  }
}