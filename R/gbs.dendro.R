#' GBS Dendrograms
#' 
#' Create different dendrograms with GBS data
#' 
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#' 
#' @param hap the hap object to manipulate
#' 
#' @keywords
#' 
#' @examples
#' 
#' @export



gbs.dendro <- function(hap) {
  
  hap.obj = hap
  
  require(ape)
  require(picante)
  require(ggplot2)
  require(rrBLUP)
  require(BLR)
  require(multicore)
  require(parallel)
  require(ape)
  
  # clean geno
  geno = read.delim(file="R/rpn/hap/rpn_geno.txt",sep="\t",header=TRUE,check.names=FALSE)
  geno = geno[grepl("SRPN",rownames(geno)),]
  
  # basic dendrogram
  hc <- hclust(dist(geno))
  
  # data frame containing labels and color divisions
  entryNames = read.delim(file="R/rpn/dendrogram/entryNames.csv",header=TRUE,check.names=FALSE,sep=",")
  entryNames = entryNames[entryNames$line_num%in%rownames(geno),]
  entryNames = (entryNames[match(hc$labels,entryNames$line_num),])
  
  # change hc labels to overlap with hc object
  summary(hc$labels==entryNames$line_num)
  
  entryNames$line_name = paste(entryNames$line_name,"-",entryNames$year)
  
  hc$labels=entryNames$line_name
  hc$labels=as.character(hc$labels)
  
  entryNames$source = factor(entryNames$source)
  entryNames$year = factor(entryNames$year)
  
  # colors
  colors = rainbow(length(levels(entryNames$source))+1,s=1,v=.9)
  colors[5] = "grey36"
  colors = colors[-(6)]
  
  # plot fan tree
  pdf("R/B.pdf",width=60,height=60)
  par(cex=1)
  op = par(bg="white")
  color.plot.phylo(as.phylo(hc), type="fan",df=entryNames,trait="source",taxa.names="line_name",leg.title="Program of Origin",label.offset=1.2,leg.cex=2,col.names=colors)
  dev.off()
  
  
  png("R/E.png",width=10,height=55,units="in",res=300)
  par(cex=.45)
  color.plot.phylo(as.phylo(hc), type="phylogram",df=entryNames,trait="year",taxa.names="line_name",label.offset=1.1,col.names=colors,no.margin = TRUE)
  dev.off()
  
  
  ######FOR STAGGERED TREE########
  # convert to dhc object
  dhc <- as.dendrogram(hc,hang=0.1)
  labels(dhc)
  
  # data frame containing labels and divisions
  
  # function to create colors based on above dataframe
  mycols <- rainbow(length(entryNames$source))[rank(entryNames$source)]
  i <- 0
  colLab <- function(n) {
    if(is.leaf(n)) {
      a <- attributes(n)
      i <<- i+1
      attr(n, "nodePar") <-
        c(a$nodePar, list(lab.col = mycols[i], lab.font = i%%3))
    }
    n
  }
  
  # apply colors to dhc object
  dL <- dendrapply(dhc, colLab)
  
  # plot staggered tree
  pdf("R/D.pdf",width=10,height=55)
  par(cex=.45)
  plot(dL, xlab="", ylab="", main="", sub="", axes=FALSE)
  dev.off()
}