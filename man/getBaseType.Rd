% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.R
\name{getBaseType}
\alias{getBaseType}
\title{Function getBaseType}
\usage{
getBaseType(dgeObj, baseType)
}
\arguments{
\item{dgeObj}{A class DGEobj created by function initDGEobj()}

\item{baseType}{One or more of: ["row", "col", "assay", "meta"]}
}
\value{
A simple list of data items
}
\description{
Accessor function for DGEobj class objects.  Retrieves all data items of a
given baseType or list of baseTypes.
}
\examples{
\dontrun{
   Assays                  <- getBaseType(DGEobj, baseType = "assay")
   AssaysAndGeneAnnotation <- getBaseType(DGEobj, c("assay", "row"))
}

}
\author{
John Thompson
}
\keyword{DGEobj}
\keyword{RNA-Seq,}
