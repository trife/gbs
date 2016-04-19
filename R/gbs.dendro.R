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