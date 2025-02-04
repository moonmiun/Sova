% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general_package_functions.r
\name{visualizeGene}
\alias{visualizeGene}
\title{Visualize a gene object with Gviz}
\usage{
visualizeGene(gene, ideogram = NA)
}
\arguments{
\item{gene}{An object of class Gene.}

\item{ideogram}{If it's different from 'NA' downloads from ensembl the
ideogram of the current chromosome. It requires internet connection and
it might require some time.}
}
\value{
A Gviz visualization of the gene object.
}
\description{
This function lets the user visualize an object of class Gene in a more
user-friendly manner.
}
\examples{
# visualizeGene(SMPD4, ideogram = "yes")

}
