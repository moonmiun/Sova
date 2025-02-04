% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes_definition.r
\docType{class}
\name{tRNA-class}
\alias{tRNA-class}
\title{Derived class: tRNA}
\description{
A class to represent tRNA
}
\section{Slots}{

\describe{
\item{\code{aminoacid}}{Aminoacid transported by the tRNA.}

\item{\code{anticodon}}{Anticodon sequence.}
}}

\examples{
# Correct notation for aminoacids
# Aminoacid information can be stored in any of the following notations:
# - full name (ex: 'phenylalanine')
# - 3 letters notation (ex: 'phe')
# - 1 letter notation (ex: 'F')
# None of these notations is required or checked by the package, but it is
# recommended for the user to only use one of them for all of the objects.
#
# Correct notation for anticodon sequences
# The standard notation for anticodon sequences is 5'-XXX-3', using RNA
# nitrogenous bases (which means Thymine is not allowed).
# For the sake of simplicity, this package only stores information about the
# basepairs and doesn't allow the user to indicate the 5' and 3' ends: it is
# good practice, though, to remember to store the anticodon information
# in the correct order.
# Anti-codons' sequences can be found here: https://rnacentral.org/.
}
