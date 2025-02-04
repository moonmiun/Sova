% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_manipulation_functions.r
\name{createGene}
\alias{createGene}
\title{Create a Gene object}
\usage{
createGene(
  id,
  symbol = NA_character_,
  name = NA_character_,
  description = NA_character_,
  assembly,
  chromosome,
  start,
  end,
  strand,
  sequence = Biostrings::DNAString("")
)
}
\arguments{
\item{id}{Gene ID.}

\item{symbol}{Gene symbol.}

\item{name}{Gene name.}

\item{description}{Description of the gene.}

\item{assembly}{Assembly used to store the position information.}

\item{chromosome}{Chromosome on which the gene is located.}

\item{start}{Start position of the gene on the chromosome (1-based).}

\item{end}{End position of the gene on the chromosome.}

\item{strand}{Strand of the gene.}

\item{sequence}{DNA sequence of the gene}
}
\value{
A Gene object.
}
\description{
This function creates an instance of the Gene class.
}
\examples{
# 0 and 1 based positions
# This package uses 1-based coordinate system since ENSEMBL uses the same
# coordinate system. The user is advised to always check the documentation of
# any database the information is collected from.
# For example, UCSC Genome Browser tables use 0-based coordinate system, which
# means that the first base pair is at position 0 and not 1. At the same time,
# the UCSC Genome Browser web interface uses the 1-start coordinate system.

# Creating a Gene object
# gene <- createGene(
#  id = "ENSG00000139618.19",
#    symbol = "BRCA2",
#    name = "BRCA2 DNA repair associated, transcript variant 4",
#   description = "Inherited mutations in BRCA1 and this gene, BRCA2, confer
#   increased lifetime risk of developing breast or ovarian cancer....",
#   assembly = "GRCh38.p14"
#   chromosome = "chr13",
#   start = 32315508,
#   end = 32400268,
#   strand = "+",
#   sequence = "ACGTACCGTCAAA....."
# )

}
