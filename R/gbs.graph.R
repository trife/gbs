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
#' 
#' @keywords
#' 
#' @examples
#' 
#' @export


gbs.graph <- function(hap,geno,file) {
  ## TODO start pdf and write all to pdf
  
  if(!missing(file)) {
    pdf(file)
  } else {
    #output to display instead
  }
  
  ## Check for blank wells, make histogram if they exist
  if(any(grepl("BLANK",colnames(hap)))) {
    missing.blank = hap[,grepl("BLANK",colnames(hap))]=="N"
    blank = as.matrix(apply(!missing.blank, 2, sum))
    
    snptot=colSums(hap[,12:ncol(hap)]!="N")
    reads=sum(snptot) #total number of reads
    sum_reads=c(snptot,blank)
    
    hist(snptot, xlim=c(0,13000), breaks=seq(0,13000, by=1000), main="SNPs number per Sample", xlab="Number of SNPs", sub="Blank wells in red")
    hist(blank, col="red",xlim=c(0,13000), breaks=seq(0,13000, by=1000),add=TRUE )
    
    hist(sum_reads/reads*100, main="% of reads for each sample", ylab="Number of Samples", xlab="Percent of total reads", sub="0.6=0.6%", xlim=c(0,0.6), breaks=seq(0, 0.6, 0.05))
    hist(blank/reads*100, add=TRUE, col="red",xlim=c(0,0.6), breaks=seq(0, 0.6, 0.05))
  }
  
  #make graphs of populations parameters
  hist(hap$MAF, main="Minor Allele Frequency", xlab="MAF Value", ylab="Number of SNPs")
  hist(hap$present, main="% Present of Each SNP", xlab="Percent Present", ylab="Number of SNPs")
  hist(hap$het, main="Number of Heterozygotes", xlab="Number of heterozygous per SNP loci", ylab="Number of SNPs")
  hist(hap$percentHET, main="Percent Heterozygous", xlab="Percent Heterozygous", ylab="Number of SNPs")
  
  #make dendrogram
  buster_dend=dendrogram(geno)
  buster_data=dendrapply(as.dendrogram(buster_dend))
  plot(buster_data, main="GBS Dendroram", type="rectangle")
  plot(as.phylo(buster_dend), type="fan", cex=0.3)
  
  #make level plot of all filtered data
  fmatch=allele.match(hap[,12:ncol(hap)])
  fmatch[lower.tri(fmatch)]=NA
  lattice::levelplot(fmatch, main="Allele matching") #TODO change color scheme

  dev.off()
}