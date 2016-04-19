#' GBS Dendrograms
#'
#' Create different dendrograms with GBS data
#'
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#' @author Trevor Rife, \email{trife@@ksu.edu}
#'
#' @param geno the geno object to manipulate
#' @param df data frame containing items to be plotted
#' @param taxa string represnting the name of the column in the dataframe that contains the names of the taxa
#' @param tips string representing the name of the column in the dataframe that contains the names for the tips
#' @param leafs string representing the name of the column in the dataframe that contains the names for the leafs
#' @param tipColor a vector of colors to be used for the tips or a column in the dataframe that should be used to divide colors
#' @param leafColor a vector of colors to be used for the leafs or a column in the dataframe that should be used to divide colors
#' @param ... additional graphical arguments from plot.phylo
#'
#' @keywords
#'
#' @examples
#'
#' @export

gbs.dendro <- function(geno, df, taxa, tips, leafs, tipColors, leafColors, ...) {

  hc <- stats::hclust(dist(geno))

  # TODO make sure hc labels and df labels match
  if(hc$labels!=df$taxa) {
    stop("Genotype labels do not match labels in input data frame.")
  }

  # TODO get colors from dataframe
  if(length(tipColor)==1) {
    colors = rainbow(length(levels(df$tipColors)))
  }

  if(length(leafColors)==1) {
    colors = rainbow(length(levels(df$tipColors)))
  }

  # TODO tips
  picante::color.plot.phylo(ape::as.phylo(hc), df=df, trait="source", taxa.names=taxa, col.names=tipColors, ...)
}

hc <- hclust(distMat)                                   # sort accessions in cluster

# load ape library for tree construction
library(ape)
hc2 <- as.phylo(hc)                                     # convert hclust object to phylo object to be used with ape

# read species info file
sppinfo <- read.table("wgrc_tetraploid_species.txt", header = T, stringsAsFactors = F)
#coreset <- as.matrix(read.table("coreset.txt", stringsAsFactors = F))

# coloring edges
edgecols=cbind(NA, rep("black", nrow(hc2$edge)), NA, NA, NA)
edgecols[,1]=hc2$tip.label[hc2$edge[,2]]                # fill in the label names using 2nd col of hc2$edge
edgecols[,1][is.na(edgecols[,1])]='N'

#edgecols[,5][edgecols[,1] %in% coreset] = "core"        # marking accessions in coreset

for (i in 1:nrow(sppinfo)){                             # filling in the spp info in 3rd column of edge object...
   for (j in 1:nrow(edgecols)){                          # based on if col 1 of sppinfo is matching col 1 of edgecols
      if (sppinfo[i,1]==edgecols[j,1]){
         edgecols[j,3]=sppinfo[i,2]
         edgecols[j,4]=sppinfo[i,3]
      }
   }
}

edgecols[,2]="black"
edgecols[,2][edgecols[,3]=="turgidum"]="red"            # coloring diff spp with different colors
edgecols[,2][edgecols[,3]=="timopheevii"]="blue"
edgecols[,2][edgecols[,4]=="durum"]="darkgreen"
edgecols[,2][edgecols[,4]=="aestivum"]="green"
#edgecols[,2][edgecols[,5]=="core"]="green"

# coloring tips
#tipcols <- cbind(hc2$tip.label)
#tipcols <- as.matrix(merge(tipcols,edgecols,by=1, sort = F))

#color minicore
mc <- read.table(file = 'minicore.txt', header = T, stringsAsFactors = F)
tipcols <- cbind(hc2$tip.label, 'black')
tipcols[,2][tipcols[,1]%in%mc$TA]='red'

# changing edge color back to only black and green
#edgecols[,2]="black"
#edgecols[,2][edgecols[,5]=="core"]="green"

library(ape)
# plotting tree
pdf("Cluster_tetraploid_w_wheat.pdf",width=42,height=42)

plot(hc2, type = 'f', lab4ut="axial", label.offset = 2, edge.color = edgecols[,2], edge.width = 3, tip.color = tipcols[,2], cex = 0.8)
title(main = substitute(paste("Cluster analysis of WGRC ", italic("Triticum turgidum"), " and ", italic("Triticum timopheevii"), " collections with wheat using GBS")), cex.main = 5, line = -3)

legend(-350,320, legend=c(expression(paste(italic("T. turgidum"))), expression(paste(italic("T. timopheevii"))), expression(paste(italic("T. turgidum ssp. durum"))), "Wheat"), text.col=c("red", "blue", "darkgreen", "green"), title.col = 'black', cex=2.5, title="Tree coloring")

legend(250,320, legend=c("Core set"), text.col=c("red"), title.col = 'black',  cex=2.5, title="Tip coloring")

dev.off()
