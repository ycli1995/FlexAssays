#' @include FlexAssays-class.R
#'
#' @importFrom methods setGeneric setMethod validObject
#' @importFrom BiocGenerics updateObject
#' @importFrom S4Vectors isTRUEorFALSE
#' @importFrom SummarizedExperiment assays assay
#' @importClassesFrom S4Vectors list_OR_List
NULL

# updateObject #################################################################

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
  assays <- updateObject(
    assays(object, withDimnames = FALSE),
    verbose = verbose
  )
  object@rowMap <- updateObject(object@rowMap, ..., verbose = verbose)
  object@colMap <- updateObject(object@colMap, ..., verbose = verbose)
  setRawAssays(object, assays)
})

# Dimensional map ##############################################################

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

# Control flexibility ##########################################################

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

# Access assays ################################################################

#' Get or set assays without dimension names
#'
#' Internal helper methods to access or update the assay list of a given
#' \code{\link{FlexAssays}}. These methods are intended for internal use during
#' development; end-users should not call these methods directly.
#'
#' @param x A \code{\link{FlexAssays}} object.
#' @param ... `r .dot_param`
#'
#' @details
#' These low-level methods are designed for internal manipulation of the assay
#' storage slot. Any subclass of \code{\link{FlexAssays}} that implements S4
#' methods for `assays` and `setRawAssays` should work well with the common
#' methods in \code{\link{FlexAssays-methods}}
#'
#' @name FlexAssays-assays
#' @keywords internal
NULL

#' @param withDimnames Logical, indicating whether to restore the dimnames of
#' assays from `rowMap` or `colMap`. Default is `TRUE`.
#'
#' @returns
#' \itemize{
#' \item `assays()`: Returns a list of matrix-like assays, optionally with
#' restored dimension names.
#' }
#'
#' @examples
#' ## assays
#' fa <- exampleFlexAssays()
#' assays(fa)
#'
#' @rdname FlexAssays-assays
#' @export
setMethod("assays", "SimpleFlexAssays", function(x, withDimnames = TRUE, ...) {
    assays <- x@data
    if (!withDimnames) {
      return(assays)
    }
    for (i in seq_along(assays)) {
      assays[[i]] <- resetDimNames(
        assays[[i]],
        rownames = mappedRowNames(x@rowMap, i),
        colnames = mappedRowNames(x@colMap, i)
      )
    }
    assays
  }
)

#' @param i Index of the assay to access or modify.
#'
#' @returns
#' \itemize{
#' \item `assay`: Returns the matrix-like assay object at the `i`-th
#' element.
#' }
#'
#' @examples
#' ## assay
#' fa <- exampleFlexAssays()
#' assay(fa, 2)
#'
#' @rdname FlexAssays-assays
#' @export
setMethod(
  "assay", c("FlexAssays", "missing"),
  function(x, i, withDimnames = TRUE, ...) {
    if (0L == length(x)) {
      stop(sprintf(
        "assay(<%s>, i=\"missing\") failed:\n  No assay exists. ",
        class(x)[1]
      ))
    }
    .get_one_assay(x, i = 1L, withDimnames = withDimnames, ...)
  }
)

#' @rdname FlexAssays-assays
#' @export
setMethod(
  "assay", c("FlexAssays", "character"),
  function(x, i, withDimnames = TRUE, ...) {
    .get_one_assay(x, i = i, withDimnames = withDimnames, ...)
  }
)

#' @rdname FlexAssays-assays
#' @export
setMethod(
  "assay", c("FlexAssays", "numeric"),
  function(x, i, withDimnames = TRUE, ...) {
    .get_one_assay(x, i = i, withDimnames = withDimnames, ...)
  }
)

.get_one_assay <- function(x, i, withDimnames = TRUE, ...) {
  assay <- assays(x, withDimnames = FALSE)[[i]]
  if (!withDimnames) {
    return(assay)
  }
  resetDimNames(assay, mappedRowNames(x@rowMap, i), mappedRowNames(x@colMap, i))
}

#' @param new.assays A list of new matrix-like assays to replace the existing
#' ones in `x`.
#' @param check Logical. If TRUE (default), [validObject()] will be called to
#' validate the object after replacement. Setting it to `FALSE` will skip checks
#' for `rowFlex`, `colFlex` and `assayClasses`.
#'
#' @returns
#' \itemize{
#' \item `setRawAssays()`: Returns an updated \code{\link{FlexAssays}} object
#' with `new.assays` set into the corresponding slot.
#' }
#'
#' @rdname FlexAssays-assays
#' @export setRawAssays
setGeneric("setRawAssays", function(x, new.assays, ..., check = TRUE) {
  standardGeneric("setRawAssays")
})

#' @examples
#' ## setRawAssays
#' fa <- exampleFlexAssays()
#' assays <- assays(fa)
#' fa <- setRawAssays(fa, assays)
#'
#' @rdname FlexAssays-assays
#' @export
setMethod(
  "setRawAssays", c("SimpleFlexAssays", "list_OR_List"),
  function(x, new.assays, ..., check = TRUE) {
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

#' @param value A single matrix-like object to set as the `i`-th assay. Must
#' contain both row and column names. If `NULL`, the original assay at index `i`
#' will be removed.
#'
#' @returns
#' \itemize{
#' \item `assay<-`: Returns an updated \code{\link{FlexAssays}} with the
#' `i`-th assay replaced by `new.assay`.
#' }
#'
#' @examples
#' ## assay<-
#' fa <- exampleFlexAssays()
#' a <- assay(fa, 1)
#' assay(fa, 4) <- a
#'
#' @rdname FlexAssays-assays
#' @export
setMethod(
  "assay<-", c("FlexAssays", "missing"),
  function(x, i, withDimnames = TRUE, ..., value) {
    if (0L == length(x)) {
      stop(sprintf(
        "`assay<-`(<%s>, i=\"missing\") failed:\n  length(assays(<%s>)) is 0. ",
        class(x)[1], class(x)[1]
      ))
    }
    .set_one_assay(x, i = 1L, new.assay = value)
  }
)

#' @rdname FlexAssays-assays
#' @export
setMethod(
  "assay<-", c("FlexAssays", "numeric"),
  function(x, i, withDimnames = TRUE, ..., value) {
    .set_one_assay(x, i = i, new.assay = value)
  }
)

#' @rdname FlexAssays-assays
#' @export
setMethod(
  "assay<-", c("FlexAssays", "character"),
  function(x, i, withDimnames = TRUE, ..., value) {
    .set_one_assay(x, i = i, new.assay = value)
  }
)

.set_one_assay <- function(x, i, new.assay) {
  assays <- assays(x, withDimnames = FALSE)
  if (is.null(new.assay)) {
    x@rowMap <- removeMapCols(x@rowMap, i)
    x@colMap <- removeMapCols(x@colMap, i)
    assays[[i]] <- new.assay
    return(setRawAssays(x, new.assays = assays, check = FALSE))
  }
  .valid_assay_classes(new.assay, x@assayClasses)
  if (nrow(new.assay) > 0 & length(rownames(new.assay)) == 0) {
    stop("'new.assay' must contain row names with non-empty rows.")
  }
  if (ncol(new.assay) > 0 & length(colnames(new.assay)) == 0) {
    stop("'new.assay' must contain column names with non-empty columns.")
  }
  if (x@rowFlex == "fixed") {
    .valid_dnames_fixed(
      new.assay,
      dimfun = rownames,
      ref = rownames(x@rowMap),
      mat.name = "new.assay"
    )
  }
  if (x@rowFlex == "bounded") {
    .valid_dnames_bounded(
      new.assay,
      dimfun = rownames,
      ref = rownames(x@rowMap),
      mat.name = "new.assay"
    )
  }
  if (x@colFlex == "fixed") {
    .valid_dnames_fixed(
      new.assay,
      dimfun = colnames,
      ref = rownames(x@colMap),
      mat.name = "new.assay"
    )
  }
  if (x@colFlex == "bounded") {
    .valid_dnames_bounded(
      new.assay,
      dimfun = colnames,
      ref = rownames(x@colMap),
      mat.name = "new.assay"
    )
  }
  row.append <- x@rowFlex == "free"
  col.append <- x@colFlex == "free"
  mappedRowNames(x@rowMap, i, append = row.append) <- rownames(new.assay)
  mappedRowNames(x@colMap, i, append = col.append) <- colnames(new.assay)
  assays[[i]] <- resetDimNames(subsetMatByDimNames(
    new.assay,
    rownames = mappedRowNames(x@rowMap, i),
    colnames = mappedRowNames(x@colMap, i),
    drop = FALSE
  ))
  setRawAssays(x, new.assays = assays, check = FALSE)
}

# Others #######################################################################

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
