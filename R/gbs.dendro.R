#' Constructs dendrogram.
#'
#' Creates a custom dendrogram with the ability to color both terminal leafs and nodes based on a separate data frame.
#'
#' @author Trevor Rife, \email{trife@@ksu.edu}
#' @author Narinder Singh, \email{nss470@@ksu.edu}
#'
#' @param phylo  An object of class phylo.
#' @param df The data frame containing the items to be plotted.
#' @param taxa A string representing the column name in the df dataframe that contains the same names that exist in the geno object.
#' @param tips A string representing the column in the df dataframe that contains the factors to use to color the tips.
#' @param leafs A string representing the column in the df dataframe that contains the factors to use to color the terminal leafs.
#' @param tipColor A vector of colors to be used for the tips.
#' @param leafColor A vector of colors to be used for the leafs.
#' @param leg.title A string representing the title for the legend.
#' @param leg.cex A numerical value for the size of the legend.
#' @param leg.color The colors to use for the legend (tipColors or leafColors).
#' @param ... Additional graphical arguments to pass to the plot.phylo function.
#'
#' @details
#' This function allows for the independent coloring of both the terminal leafs of a dendrogram as well as the tips based on an a separate data frame.
#'
#' @examples
#' data(wheat)
#' phylo <- as.phylo(hclust(dist(hap$geno)))
#' gbs.dendro(phylo,phenotypes,taxa="line",tips="type",leafs="group",type="fan",edge.width=2)
#'
#' @export

gbs.dendro <- function(phylo, df, taxa, tips, leafs, tipColors, leafColors, leg.title = NULL, leg.cex = 1, leg.color = tipColors, ...) {
  
  # Geno data integrity
  if (class(phylo) != "phylo") {
    stop("phylo must be a phylo object")
  }
  
  hc2 <- phylo
  
  # Data integrity
  len.tips <- length(hc2$tip.label)
  len.taxa <- length(df[[taxa]])
  
  if (len.tips != len.taxa || sum(hc2$tip.label %in% df[[taxa]]) != len.taxa) {
    stop("Missing taxa in tree or data frame")
  }
  
  # Color tips
  tip.color <- cbind(hc2$tip.label, 'black')
  
  if(missing(tips)) {
    tip.color <- rep("black",length(hc2$tip.label))
  } else {
    if(missing(tipColors)) {
      tipColors <- rainbow(length(levels(df[[tips]])))
    }
    
    # Match order of the data frame matches the tips
    order <- match(hc2$tip.label, df[,taxa])
    ordered.trait <- df[tips][order,]
    
    # Cut up the tip trait and assign a list of colors
    levs <- levels(ordered.trait)
    tip.color <- rep("black", times = length(df[[taxa]]))
    tip.color <- tipColors[match(ordered.trait, levs)]
  }
  
  if(missing(leafs)) {
    edgeCols <- cbind(NA, rep("black", nrow(hc2$edge)), NA, NA)
  } else {
    edgeCols <- cbind(NA, rep("black", nrow(hc2$edge)), NA, NA)
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
      leafColors <- rainbow(length(levels(df[[leafs]])))
    }
    
    # Match the order of the data frame matches the leafs
    order <- match(edgeCols[edgeCols[,1]!="N",][,1], df[,taxa])
    ordered.trait <- df[leafs][order,]
    
    # Cut up the leaf trait and assign a list of colors
    levs <- levels(ordered.trait)
    edge.color <- rep("black", times = length(df[[taxa]]))
    edge.color <- leafColors[match(ordered.trait, levs)]
    
    edgeCols[edgeCols[,1]!="N",][,2]=edge.color
  }

  # Plot
  ape::plot.phylo(hc2, tip.color = tip.color, edge.color = edgeCols[,2], ...)
  
  legend("bottomleft", levs, fill = leg.color, inset = 0.05, title = leg.title, cex = leg.cex)
}