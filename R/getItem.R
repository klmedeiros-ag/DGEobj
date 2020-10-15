#' Function getItems
#'
#' Retrieve items from a DGEobj by item name.
#'
#' @author John Thompson
#' @keywords RNA-Seq, DGEobj
#'
#' @param dgeObj  A class DGEobj created by function initDGEobj()
#' @param itemNames A list of itemNames to retrieve
#'
#' @return The requested data item or list of data items.
#'
#' @examples
#' \dontrun{
#'    MyCounts <- getItem(DGEobj, "counts")
#' }
#'
#' @importFrom assertthat assert_that
#'
#' @export
getItems <- function(dgeObj, itemNames){

    assertthat::assert_that(!missing(dgeObj),
                            !missing(itemNames),
                            msg = "Be sure to specify a DGEobj and at least one itemName to retrieve.")
    assertthat::assert_that("DGEobj" %in% class(dgeObj),
                            msg = "The DGEobj must be of class 'DGEobj'.")
    assertthat::assert_that(any(c("character", "list") %in% class(itemNames)),
                            msg = "Pass the itemNames as a single character string or a list of items to retrieve.")

    idx <- itemNames %in% names(dgeObj)
    result <- list()
    for (itemName in itemNames[idx]) {
        result[[itemName]] <- dgeObj[[itemName]]
    }

    if (length(result) == 1) result <- result[[1]]

    if (sum(idx) < length(idx)) {
        missingItems <- stringr::str_c(itemNames[!idx], sep = ", ")
        warning(stringr::str_c("These item(s) not found: [", missingItems, "]"))
    }

    return(result)
}


#' Function getItem
#'
#' Retrieve an item from a DGEobj by item name.
#'
#' @author John Thompson
#' @keywords RNA-Seq, DGEobj
#'
#' @param dgeObj  A class DGEobj created by function initDGEobj()
#' @param itemName Name of item to retrieve
#'
#' @return The requested data item
#'
#' @examples
#' \dontrun{
#'    MyCounts <- getItem(DGEobj, "counts")
#' }
#'
#' @importFrom assertthat assert_that
#'
#' @export
getItem <- function(dgeObj, itemName){
    assertthat::assert_that(!missing(dgeObj),
                            !missing(itemName),
                            msg = "Be sure to specify a DGEobj and an itemName to retrieve.")
    assertthat::assert_that("DGEobj" %in% class(dgeObj),
                            msg = "The DGEobj must be of class 'DGEobj'.")
    assertthat::assert_that(class(itemName) == "character",
                            length(itemName) == 1,
                            msg = "The itemName should be a character string and contain the name of only one item to retrieve.
                             To retrieve multiple items, use the getItems() function.")
    assertthat::assert_that(itemName %in% names(dgeObj),
                            msg = "The requested itemName should be in the DGEobj. Use names(dgeObj) to see the available items.")
    return(dgeObj[[itemName]])
}


#' Function getType
#'
#' Retrieve one or more data items from a DGEobj by type.
#'
#' @author John Thompson
#' @keywords RNA-Seq, DGEobj
#'
#' @param dgeObj  A class DGEobj created by function initDGEobj()
#' @param type A single type of list of types to retrieve.  Enter
#'    showTypes(MyDGEobj) to see a list of allowed types.  See addType() function
#'    to define new types.
#' @param parent (Optional) Filter return list for common parent (e.g. useful
#' to select one set of contrast results when multiple fits have been performed)
#'
#' @return A list of requested data items
#'
#' @examples
#' \dontrun{
#'    MyContrastList <- getType(DGEobj, type = "topTable")
#'    MyRawData      <- getType(DGEobj, type = list("counts", "design", "geneData"))
#'}
#'
#' @export
getType <- function(dgeObj, type, parent){

    idx <- attr(dgeObj, "type") %in% type
    if (!missing(parent)) {
        pidx <- attr(dgeObj, "parent") == parent
        idx <- idx & pidx
    }
    result <- unclass(dgeObj)[idx]

    if (sum(idx) < length(type))
        warning("Some types were not found")
    if (sum(idx) == 0) {
        tsmsg("Warning: no items of specified type are found.")
        return(NULL)
    } else {
        if (sum(idx) < length(type))
            warning("Some types were not found")
        return(result)
    }
}


#' Function getBaseType
#'
#' Accessor function for DGEobj class objects.  Retrieves all data items of a
#' given baseType or list of baseTypes.
#'
#' @author John Thompson
#' @keywords RNA-Seq, DGEobj
#'
#' @param dgeObj  A class DGEobj created by function initDGEobj()
#' @param baseType One or more of: ["row", "col", "assay", "meta"]
#'
#' @return A simple list of data items
#'
#' @examples
#' \dontrun{
#'    Assays                  <- getBaseType(DGEobj, baseType = "assay")
#'    AssaysAndGeneAnnotation <- getBaseType(DGEobj, c("assay", "row"))
#' }
#'
#' @export
getBaseType <- function(dgeObj, baseType){

    if (missing(baseType))
        stop("baseType argument is required")

    if (!baseType %in% baseTypes(dgeObj))
        stop(paste("baseType must be one of: ",
                   paste(baseTypes(dgeObj), collapse = ", "),
                   sep = ""))

    idx <- attr(dgeObj, "basetype") %in% baseType

    if (sum(idx) < length(baseType))
        warning("Some baseTypes were not found")

    result <- unclass(dgeObj)[idx]
    return(result)
}


#' Function BaseType
#'
#' Return the baseType for a given item type in a DGEobj.
#'
#' @author John Thompson
#' @keywords RNA-Seq, DGEobj
#'
#' @param dgeObj A class DGEobj created by function initDGEobj()
#' @param type  An item type for which you want the baseType
#'
#' @return A baseType value (character string)
#'
#' @examples
#' \dontrun{
#'    MyBaseType <- baseType(dgeObj, type = "DGEList")
#' }
#'
#' @importFrom assertthat assert_that
#'
#' @export
baseType <- function(dgeObj, type){

    assertthat::assert_that(!missing(dgeObj),
                            !missing(type),
                            msg = "Be sure to specify a DGEobj and an item type.")
    assertthat::assert_that(class(dgeObj)[[1]] == "DGEobj",
                            msg = "The DGEobj must be of class 'DGEobj'.")
    assertthat::assert_that(class(type)[[1]] == "character",
                            msg = "The type must be of class 'character'.")

    objDef <- attr(dgeObj, "objDef")
    return(objDef$type[[type]])
}


#' Function baseTypes
#'
#' Return a list of the available baseTypes in a DGEobj.
#'
#' @author John Thompson
#' @keywords RNA-Seq, DGEobj
#'
#' @param dgeObj  A class DGEobj object
#'
#' @return A list of baseTypes
#'
#' @examples
#' \dontrun{
#'    # Global definition of baseTypes
#'    myBaseTypes <- baseTypes()
#'
#'    # Basetypes from a specific DGEobj
#'    myBaseTypes <- baseTypes(myDGEobj)
#' }
#'
#' @export
baseTypes <- function(dgeObj){
    if (missing(dgeObj))
        return(unique(.DGEobjDef$type))
    else
        return(unique(attr(dgeObj, "objDef")$type))
}
