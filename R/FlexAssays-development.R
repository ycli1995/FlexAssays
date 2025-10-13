#' @include utils.R
#' @include logical_map.R
#' @include FlexAssays.R
#'
#' @importFrom methods setGeneric setMethod validObject
#' @importFrom BiocGenerics updateObject
#' @importFrom S4Vectors isTRUEorFALSE
#' @importFrom SummarizedExperiment assay assay<- assays assays<-
#' @importClassesFrom S4Vectors list_OR_List
NULL

#' Update a FlexAssays object
#'
#' Method to update \code{\link{FlexAssays}} objects to the latest version.
#'
#' @param object A legacy \code{\link{FlexAssays}} object that may need to be
#' updated.
#' @param verbose Logical. If `TRUE`, messages about the update process will be
#' printed.
#' @param ... `r .dot_param`
#'
#' @return
#' An updated version of `object` compatible with the current `FlexAssays` class
#' structure.
#'
#' @examples
#' fa <- exampleFlexAssays()
#' updateObject(fa)
#'
#' @export
#' @aliases updateObject updateObject,FlexAssays-method
#' @name FlexAssays-updateObject
setMethod("updateObject", "FlexAssays", function(object, ..., verbose = FALSE) {
  a <- updateObject(assays(object, withDimnames = FALSE), verbose = verbose)
  object@rowMap <- updateObject(object@rowMap, ..., verbose = verbose)
  object@colMap <- updateObject(object@colMap, ..., verbose = verbose)
  setRawAssays(object, a, check = FALSE)
})

#' Dimensional map for FlexAssays
#'
#' Accessors for the row/column mapping matrices of a \code{\link{FlexAssays}}
#' object. These logical maps indicate which rows or columns are present in each
#' assay. Primarily intended for developers rather than end-users.
#'
#' @param x A \code{\link{FlexAssays}} object
#'
#' @returns
#' A logical `r .doc_links("SVT_SparseMatrix")` indicating the presence of
#' global row or column names across individual assays in `x`.
#'
#' @examples
#' fa <- exampleFlexAssays()
#' rowMap(fa)
#' colMap(fa)
#'
#' @name FlexAssays-map
NULL

#' @export
#' @rdname FlexAssays-map
setGeneric("rowMap", function(x) standardGeneric("rowMap"))

#' @export
#' @rdname FlexAssays-map
setMethod("rowMap", "FlexAssays", function(x) x@rowMap)

#' @export
#' @rdname FlexAssays-map
setGeneric("colMap", function(x) standardGeneric("colMap"))

#' @export
#' @rdname FlexAssays-map
setMethod("colMap", "FlexAssays", function(x) x@colMap)

#' Control Flexibility for FlexAssays
#'
#' Methods for retrieving the configuration parameters that control row/column
#' flexibility and supported assay classes in a \code{\link{FlexAssays}} object.
#'
#' @param x A \code{\link{FlexAssays}} object
#'
#' @details
#' The flexibility of rows and columns in \code{\link{FlexAssays}} is governed
#' by one of the following rules:
#' \itemize{
#' \item `free`: Global row/column names can be extended when adding assays with
#' new dimension names.
#' \item `bounded`: Global row/column names are fixed; assays with unmatched
#' names ***cannot*** be added.
#' \item `fixed`: All assays must share the exact same row/column names.
#' }
#'
#' @examples
#' fa <- exampleFlexAssays()
#' rowFlex(fa)
#' colFlex(fa)
#' assayClasses(fa)
#'
#' @name FlexAssays-flex
NULL

#' @returns
#' \itemize{
#' \item `rowFlex`, `colFlex`: A character scalar indicating the rule of row or
#' column flexibility.
#' }
#'
#' @export
#' @rdname FlexAssays-flex
setGeneric("rowFlex", function(x) standardGeneric("rowFlex"))

#' @export
#' @rdname FlexAssays-flex
setMethod("rowFlex", "FlexAssays", function(x) x@rowFlex)

#' @export
#' @rdname FlexAssays-flex
setGeneric("colFlex", function(x) standardGeneric("colFlex"))

#' @export
#' @rdname FlexAssays-flex
setMethod("colFlex", "FlexAssays", function(x) x@colFlex)

#' @returns
#' \itemize{
#' \item `assayClasses`: A character vector indicating the supported classes for
#' assay elements.
#' }
#'
#' @export
#' @rdname FlexAssays-flex
setGeneric("assayClasses", function(x) standardGeneric("assayClasses"))

#' @export
#' @rdname FlexAssays-flex
setMethod("assayClasses", "FlexAssays", function(x) x@assayClasses)

#' Update all assays with/without validation
#'
#' These functions are intended for low-level manipulation, e.g.,
#' programmatically replacing all matrix-like objects in a sub-class of
#' [`FlexAssays`].
#'
#' @name setRawAssays
NULL

#' @param x A [`FlexAssays`] object.
#' @param ... `r .dot_param`
#'
#' @details
#' `setRawAssays` is a lower-level function than [`assays<-`] and is only
#' intended for development usage. It can bypass the time consuming validation,
#' therefore the caller must ensure `new.assays` are fitted with the correct
#' slot for the subclass of [`FlexAssays`].
#'
#' @returns
#' `setRawAssays` returns an updated [`FlexAssays`] object with `new.assays` set
#' into the corresponding slot.
#'
#' @rdname setRawAssays
#' @export
setGeneric("setRawAssays", function(x, new.assays, check = TRUE, ...) {
  standardGeneric("setRawAssays")
})

#' @param new.assays A list of new matrix-like assays to replace the existing
#' ones in `x`. If is `NULL`, should remove all internal assays.
#' @param check Logical. If `TRUE` (default), [validObject()] will be called to
#' validate the object after replacement. Setting it to `FALSE` will skip checks
#' for `rowFlex`, `colFlex` and `assayClasses`.
#'
#' @examples
#' fa <- exampleFlexAssays()
#' aa <- assays(fa)[1:2]
#'
#' ## ERROR: The number of assays is 2, but the dimensional map has 3 columns.
#' try(fa <- setRawAssays(fa, aa))
#'
#' @rdname setRawAssays
#' @export
setMethod(
  "setRawAssays", c("SimpleFlexAssays", "list_OR_List"),
  function(x, new.assays, check = TRUE, ...) {
    stopifnot(isTRUEorFALSE(check))
    for (i in seq_along(new.assays)) {
      new.assays[[i]] <- resetDimNames(new.assays[[i]])
    }
    x@data <- as(new.assays, "SimpleList")
    if (check) {
      validObject(x)
    }
    x
  }
)

#' @returns
#' `cleanRawAssays` returns an updated [`FlexAssays`] object with the internal
#' assay data and the corresponding dimensional maps are cleaned.
#'
#' @rdname setRawAssays
#' @export
setGeneric("cleanRawAssays", function(x, check = TRUE, ...) {
  standardGeneric("cleanRawAssays")
})

#' @examples
#' ## Clean all assays while keep the dimensional rules.
#' fa <- exampleFlexAssays()
#' fa <- cleanRawAssays(fa)
#'
#' @rdname setRawAssays
#' @export
setMethod("cleanRawAssays", "SimpleFlexAssays", function(x, check = TRUE, ...) {
  x@rowMap <- x@rowMap[, integer(), drop = FALSE]
  x@colMap <- x@colMap[, integer(), drop = FALSE]
  x@data <- SimpleList()
  if (check) {
    validObject(x)
  }
  x
})

#' Clean Unmapped Dimensions for FlexAssays
#'
#' Removes unmapped row and column names from the `rowMap` and `colMap` of a
#' \code{\link{FlexAssays}} object, effectively trimming unused dimensions. This
#' function is primarily intended for developers or pipeline maintainers to
#' optimize memory and ensure consistency in dimension maps.
#'
#' @param x A \code{\link{FlexAssays}} object
#' @param check Logical. If `TRUE` (default), [validObject()] will be called
#' after cleaning to ensure the object remains valid.
#' @param ... `r .dot_param`
#'
#' @returns
#' A \code{\link{FlexAssays}} object with cleaned dimension maps.
#'
#' @examples
#' fa <- exampleFlexAssays()
#' fa[[2]] <- NULL
#'
#' fa <- cleanDimMaps(fa)
#'
#' @export
cleanDimMaps <- function(x, check = TRUE, ...) {
  stopifnot(inherits(x, "FlexAssays"))
  x@rowMap <- dropMapRows(x@rowMap)
  x@colMap <- dropMapRows(x@colMap)
  if (check) {
    validObject(x)
  }
  x
}
