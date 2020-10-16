#' Function annotateDGEobj
#'
#' Reads an annotation file  of key/value pairs and attach them attributes to a DGEobj.
#' The annotation file should be a text file containing key/value pairs separated by an equals
#' sign.  The keys parameter specifies which key we want to capture as
#' attributes on the DGEobj. The value will then be the value of that attribute.
#'
#' @author John Thompson
#' @keywords RNA-Seq, DGEobj, annotation, attributes
#'
#' @param dgeObj  A class DGEobj created by function initDGEobj()
#' @param annotationFile An Omicsoft registration file (in text format, not Excel)
#' @param keys A list of keys to look for in the annotationFile and transfer to the DGEobj.
#'    Use keys = NULL to accept all keys in the annotationFile (Default)
#'
#' @return A DGEobj annotated with attributes from the annotation file.
#'
#' @examples
#' \dontrun{
#'    MyDgeObj <- annotateDGEobj(DGEobj, annotationFile)
#' }
#'
#' @import magrittr
#' @importFrom stringr str_detect str_remove_all str_locate
#' @importFrom utils read.delim
#'
#' @export
annotateDGEobj <- function(dgeObj, annotationFile, keys = NULL) {

    assertthat::assert_that(file.exists(annotationFile),
                            msg = "You must provide an annotation text file (annotationFile) which contains key/value pairs separated by an equals sign.")

    # Read lines, stripping quotes
    regdat <- utils::read.delim(annotationFile, sep = "\t",
                                quote = "\"",
                                stringsAsFactors = FALSE,
                                header = FALSE)
    # Just first column
    regdat <- regdat[,1, drop = FALSE]
    colnames(regdat) <- "pair"

    # Just lines with equals signs
    regdat <- regdat[grepl("=", regdat$pair) | !grepl("Parameters.", regdat$pair), , drop = FALSE]

    # Loop through the attributes spitting on the first = sign
    regdat$key <- ""
    regdat$value <- ""
    for (i in 1:nrow(regdat)) {
        splitpos <- str_locate(regdat[i,1], "=")[1]  # Pos of 1st = sign
        regdat$key[i] <- substr(regdat[i,1], 1, (splitpos - 1))
        regdat$value[i] <- substr(regdat[i,1], (splitpos + 1), nchar(regdat[i,1]))
    }

    # After splitting, key without values get the key names inserted as the value. Convert those to empty strings.
    idx <- regdat$key == regdat$value
    regdat$value[idx] <- ""

    # Squeeze spaces out of keys
    regdat$key <- stringr::str_remove_all(regdat$key, " ")

    # Capture/preserve the existing attributes
    MyAttribs <- attributes(dgeObj)
    if (is.null(MyAttribs)) MyAttribs <- list()

    for (i in 1:nrow(regdat))
        if (is.null(keys)) { # Add all keys
            MyAttribs[regdat$key[i]] <- regdat$value[i]
        } else {
            # Just add specified keys
            if (regdat$key[i] %in% keys)
                MyAttribs[regdat$key[i]] <- regdat$value[i]
        }

    # Now attach the attributes to the DGEobj
    attributes(dgeObj) <- MyAttribs
    return(dgeObj)
}
