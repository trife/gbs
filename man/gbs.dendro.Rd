% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbs.dendro.R
\name{gbs.dendro}
\alias{gbs.dendro}
\title{Constructs dendrogram.}
\usage{
gbs.dendro(phylo, df, taxa, tips, leafs, tipColors, leafColors,
  leg.title = NULL, leg.cex = 1, leg.color = tipColors, ...)
}
\arguments{
\item{phylo}{An object of class phylo.}

\item{df}{The data frame containing the items to be plotted.}

\item{taxa}{A string representing the column name in the df dataframe that contains the same names that exist in the geno object.}

\item{tips}{A string representing the column in the df dataframe that contains the factors to use to color the tips.}

\item{leafs}{A string representing the column in the df dataframe that contains the factors to use to color the terminal leafs.}

\item{leg.title}{A string representing the title for the legend.}

\item{leg.cex}{A numerical value for the size of the legend.}

\item{leg.color}{The colors to use for the legend (tipColors or leafColors).}

\item{...}{Additional graphical arguments to pass to the plot.phylo function.}

\item{tipColor}{A vector of colors to be used for the tips.}

\item{leafColor}{A vector of colors to be used for the leafs.}
}
\description{
Creates a custom dendrogram with the ability to color both terminal leafs and nodes based on a separate data frame.
}
\details{
This function allows for the independent coloring of both the terminal leafs of a dendrogram as well as the tips based on an a separate data frame.
}
\examples{
data(wheat)
phylo <- as.phylo(hclust(dist(hap$geno)))
gbs.dendro(phylo,phenotypes,taxa="line",tips="type",leafs="group",type="fan",edge.width=2)

}
\author{
Trevor Rife, \email{trife@ksu.edu}

Narinder Singh, \email{nss470@ksu.edu}
}

