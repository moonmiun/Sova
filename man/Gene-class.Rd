% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes_definition.r
\docType{class}
\name{Gene-class}
\alias{Gene-class}
\title{Class Gene}
\description{
A virtual class representing a generic gene.
}
\section{Slots}{

\describe{
\item{\code{id}}{A unique identifier for the gene (ex, Ensembl ID 'ENSG00000139618').}

\item{\code{symbol}}{Gene symbol (ex, HUGO symbol like 'BRCA2').}

\item{\code{name}}{Complete name for the gene.}

\item{\code{description}}{Gene description. Often the summary found on NCBI database.}

\item{\code{assembly}}{Assembly used to store the position information (better to use
ENSEMBL IDs 'hg19' or 'hg38').}

\item{\code{chromosome}}{Chromosome position.}

\item{\code{start}}{Starting position (0 based).}

\item{\code{end}}{Ending position.}

\item{\code{strand}}{DNA strand ('+' or '-').}

\item{\code{sequence}}{DNA sequence}
}}

\keyword{classes}
