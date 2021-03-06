% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbs.summary.R
\name{gbs.summary}
\alias{gbs.summary}
\title{summarize the gbs object.}
\usage{
gbs.summary(hap, geno = F, encoding = c(-1, 0, 1), graph = F, het = "H",
  missing = "N")
}
\arguments{
\item{hap}{The gbs object to be processed.}

\item{geno}{A logical value that will convert the marker calls to a numeric geno format.}

\item{encoding}{The numbers used to encode alleles (minor, het, major).}

\item{graph}{A logical value that will output graphs for summary statistics (blank wells, maf, percent present, het, percent het).}

\item{het}{The symbol(s) used for heterozygous calls.}

\item{missing}{The symbol(s) used for missing calls.}
}
\description{
Processes the gbs object and calculates allele statistics. Returns allele stats, a geno object, and graphs statistics.
}
\details{
This function checks for and removes any wells that contain "blank" in the sample name, calculates allele counts, percent present, and minor allele frequency, and optionally creates a numeric geno object (stored as hap$geno). The encoding argument can be used to change how the minor, major, and heterozygous alleles are encoded. The het and missing arguments can be a character or vector to allow for multiple matches (e.g. IUPAC het = c("R","Y","S","W","K","M")).
}
\examples{
data(wheat)
hap = gbs.summary(hap,geno=T,graph=T)

}
\author{
Trevor Rife, \email{trife@ksu.edu}

Jesse Poland, \email{jpoland@ksu.edu}
}

