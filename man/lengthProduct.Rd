% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general_package_functions.r
\name{lengthProduct}
\alias{lengthProduct}
\alias{lengthProduct,Gene-method}
\alias{lengthProduct,ProteinCodingGene-method}
\alias{lengthProduct,ncRNA-method}
\alias{lengthProduct,miRNA-method}
\title{Length of the genetic product function}
\usage{
lengthProduct(object)

\S4method{lengthProduct}{Gene}(object)

\S4method{lengthProduct}{ProteinCodingGene}(object)

\S4method{lengthProduct}{ncRNA}(object)

\S4method{lengthProduct}{miRNA}(object)
}
\arguments{
\item{object}{Object of class Gene and its derivates.}
}
\value{
The length of the genetic product.
}
\description{
Calculates the length of the genetic product of each class. What is the
genetic product differs from class to class.
}
\examples{
# lengthProduct(SMPD4)

# Class Object
# This is the case in which there is no information about the gene so the
# genetic product is the length of its RNA transcription.
# lengthProduct(Gene)
# Class ProteinCodingGene
# Example in which the RNA sequence is given.
# set_slot(protein_gene, 'protein_sequence', 'MTEYKLVVVG')
# lengthProduct(protein_gene)
# > 10

# Example in which the RNA sequence is not given so the function uses the
# length of the exons.
# lengthProduct(protein_gene)
# Class ncRNA
# In this case the genetic product is simply the RNA transcription of the gene,
# which has the same length as the gene itself.
# This method is valid also for objects of sub-classes
# 'tRNA', 'rRNA' and 'snRNA'.
# lengthProduct(ncRNA)
# Class miRNA
# set_slot(miRNA_object, 'mature_sequence', 'UGAGGUAGUAGGUUGUAUAGUU')
# lengthProduct(miRNA_object)
# > 22
}
