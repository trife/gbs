#' GBS Dendrograms
#'
#' Create different dendrograms with GBS data borrows some code from picante package
#'
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#'
#' @param geno the geno object to manipulate
#' @param df data frame containing items to be plotted
#' @param taxa string representing the column in the dataframe that contains the names of the taxa in the geno object
#' @param tips string representing the column in the dataframe that contains the trait of interest for the tips
#' @param leafs string representing the column in the dataframe that contains the trait of interest for the leafs
#' @param tipColor a vector of colors to be used for the tips
#' @param leafColor a vector of colors to be used for the leafs
#' @param ... additional graphical arguments to pass to the plot.phylo function
#'
#' @keywords
#'
#' @examples
#'
#' @export

gbs.dendro <- function(geno, df, taxa, tips, leafs, tipColors, leafColors, ...) {
  
  # Convert hclust object to phylo object
  hc <- stats::hclust(dist(geno))
  hc2 <- ape::as.phylo(hc)
  
  # Data integrity
  len.tips <- length(hc2$tip.label)
  len.taxa <- length(df[[taxa]])
  
  if (len.tips != len.taxa | sum(hc2$tip.label %in% df[[taxa]]) != len.taxa) {
    stop("Missing taxa in tree or data frame")
  }
  
  # color tips
  tip.color <- cbind(hc2$tip.label, 'black')
  
  if(missing(tipColors)) {
    tip.color = rainbow(length(levels(df[[tips]])))
  }
  
  # match order of the data frame matches the tips
  order <- match(hc2$tip.label, df[,taxa])
  ordered.trait <- df[tips][order,]
  
  # cut up the tip trait and assign a list of colors
  levs <- levels(ordered.trait)
  tip.color <- rep("black", times = length(df[[taxa]]))
  tip.color <- tipColors[match(ordered.trait, levs)]
  
  
  # color edges
  edgeCols=cbind(NA, rep("black", nrow(hc2$edge)), NA, NA)
  edgeCols[,1]=hc2$tip.label[hc2$edge[,2]]
  edgeCols[,1][is.na(edgeCols[,1])]='N'
  
  for (i in 1:nrow(df)){
    for (j in 1:nrow(edgeCols)){
      if (df[i,1]==edgeCols[j,1]){
        edgeCols[j,3]=as.character(df[i,2])
        edgeCols[j,4]=as.character(df[i,3])
      }
    }
  }
  
  if(missing(leafColors)) {
    leafColors = rainbow(length(levels(df[[leafs]])))
  }
  
  # match the order of the data frame matches the leafs
  order <- match(edgeCols[edgeCols[,1]!="N",][,1], df[,taxa])
  ordered.trait <- df[leafs][order,]
  
  # cut up the leaf trait and assign a list of colors
  levs <- levels(ordered.trait)
  edge.color <- rep("black", times = length(df[[taxa]]))
  edge.color <- leafColors[match(ordered.trait, levs)]
  
  edgeCols[edgeCols[,1]!="N",][,2]=edge.color
  
  # plot
  ape::plot.phylo(hc2, tip.color = tip.color, edge.color = edgeCols[,2], ...)
}