% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addItem.R
\name{addItem}
\alias{addItem}
\title{Function addItem (DGEobj)}
\usage{
addItem(
  dgeObj,
  item,
  itemName,
  itemType,
  overwrite = FALSE,
  funArgs = match.call(),
  itemAttr,
  parent = "",
  init = FALSE,
  debug = FALSE
)
}
\arguments{
\item{dgeObj}{A class DGEobj created by function initDGEobj()}

\item{item}{The data item to be deposited in the DGEobj (Required)}

\item{itemName}{A user assigned name for this data item. (Required)}

\item{itemType}{A type attribute.  See showTypes() to see the
predefined types. Types are extensible with the newType() function. (Required)}

\item{overwrite}{Default = FALSE.  Set to TRUE to overwrite the data object
stored in the itemName slot}

\item{funArgs}{A text field to annotate how the data object was created.
If you pass the result of match.call() as this argument, it captures the
name and arguments used in the current function (optional)}

\item{itemAttr}{A named list of attributes to add directly to the item (optional)}

\item{parent}{itemName of the parent of this item (optional, but your DGEobj
won't be well annotated if you don't use this wherever appropriate)}

\item{init}{Default = FALSE. Used internally by the initDGEobj() function.}

\item{debug}{Default = FALSE; TRUE trigger browser mode.}
}
\value{
A DGEobj class object with a new data item added.
}
\description{
Add a data item to a class DGEobj
}
\examples{
\dontrun{
   myFunArgs <- match.call() # Capture calling function and arguments
   showTypes()  # See what predefined types are available
   myDGEobj <- addItem(myDGEobj, item = MyCounts,
                                 itemName = "counts",
                                 itemType = "counts",
                                 funArgs = myFunArgs)
}

}
\author{
John Thompson
}
\keyword{DGEobj}
\keyword{RNA-Seq,}
