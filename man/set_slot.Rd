% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_manipulation_functions.r
\name{set_slot}
\alias{set_slot}
\title{Set a specific slot value}
\usage{
set_slot(object, slot_name, value)
}
\arguments{
\item{object}{Object of Gene, and its derivates, class.}

\item{slot_name}{Name of the slot.}

\item{value}{The new value to assign to the slot.}
}
\value{
The updated object.
}
\description{
Set a slot value without directly having to access it; the value can be
either added or edited.
}
\examples{
# Changing phenylalanine notation for the object tRNA created before.
# tRNA <- set_slot(tRNA, 'aminoacid', 'Phe')

}
