#' @include FlexAssays-methods.R
#'
#' @importFrom utils .DollarNames
NULL

#' Subset a FlexAssays
#'
#' Provides methods for extracting or subsetting elements from a
#' \code{\link{FlexAssays}} object with indices or names.
#'
#' @param x A \code{\link{FlexAssays}} object.
#' @param ... `r .dot_param`
#'
#' @seealso [FlexAssays()] for object construction and structure.
#'
#' @name FlexAssays-subset
NULL

#' @param i,j Indices specifying elements to extract or replace. In `[[` or
#' `[[<-`, argument `j` is ignored.
#'
#' @returns
#' \itemize{
#' \item `[[`: Returns a matrix-like assay with dimension names for the `i`-th
#' element of `x`.
#' }
#'
#' @export
#' @rdname FlexAssays-subset
setMethod(
  f = "[[",
  signature = c("FlexAssays", "ANY", "missing"),
  definition = function(x, i, j, ...) assay(x, i)
)

#' @param value `r .val_param` Typically a matrix-like object to be set. If is
#' `NULL`, the corresponding assay will be removed.
#'
#' @returns
#' \itemize{
#' \item `[[<-`: Returns an updated `FlexAssays` object with the `i`-th assay
#' modified.
#' }
#'
#' @examples
#' fa <- exampleFlexAssays()
#' fa[[4]] <- fa[[1]]
#'
#' @aliases [[<-,FlexAssays,ANY,missing-method
#'
#' @export
#' @rdname FlexAssays-subset
setMethod(
  f = "[[<-",
  signature = c("FlexAssays", "ANY", "missing"),
  definition = function(x, i, j, ..., value) {
    assay(x, i = i) <- value
    x
  }
)

#' @param pattern A regular expression used in `.DollarNames` to filter names.
#' Only matching names are returned.
#'
#' @returns
#' \itemize{
#' \item `.DollarNames`: Returns a character vector of assay names matching
#' `pattern`, used for tab-completion.
#' }
#'
#' @seealso [.DollarNames()] for the S3 generic.
#'
#' @method .DollarNames FlexAssays
#' @export
#' @rdname FlexAssays-subset
.DollarNames.FlexAssays <- function(x, pattern = "") {
  grep(pattern, names(assays(x, withDimnames = FALSE)), value = TRUE)
}

#' @param name The name of the assay to get or set using the `$` accessor.
#'
#' @returns
#' \itemize{
#' \item `$`: Returns a matrix-like assay with the specified `name` in `x`.
#' }
#'
#' @rdname FlexAssays-subset
#' @export
#' @method $ FlexAssays
"$.FlexAssays" <- function(x, name) assay(x, name)

#' @returns
#' \itemize{
#' \item `$<-`: Returns an updated `FlexAssays` object with a new or modified
#' assay.
#' }
#'
#' @rdname FlexAssays-subset
#' @export
#' @method $<- FlexAssays
"$<-.FlexAssays" <- function(x, name, value) {
  assay(x, i = name) <- value
  x
}

#' @param drop Logical, whether or not to drop those empty layers after
#' subsetting.
#'
#' @returns
#' \itemize{
#' \item `[`: Returns a subset of the original `FlexAssays` object, optionally
#' dropping empty assays if `drop = TRUE`.
#' }
#'
#' @export
#' @rdname FlexAssays-subset
setMethod("[", "FlexAssays", function(x, i, j, ..., drop = TRUE) {
  miss.i <- missing(i)
  miss.j <- missing(j)
  if (miss.i & miss.j) {
    return(x)
  }
  rowmap <- x@rowMap
  colmap <- x@colMap
  if (!miss.i) {
    rowmap <- rowmap[i, , drop = FALSE]
  }
  if (!miss.j) {
    colmap <- colmap[j, , drop = FALSE]
  }
  if (any(nrow(rowmap) == 0, nrow(colmap) == 0) & drop) {
    x@rowMap <- rowmap[, integer(), drop = FALSE]
    x@colMap <- colmap[, integer(), drop = FALSE]
    return(setRawAssays(x, NULL, check = FALSE))
  }
  keep <- logical(length(x))
  assays <- assays(x)
  for (i in seq_along(assays)) {
    assays[[i]] <- resetDimNames(subsetMatByDimNames(
      assays[[i]],
      mappedRowNames(rowmap, i),
      mappedRowNames(colmap, i),
      drop = FALSE
    ))
    if (prod(dim(assays[[i]])) > 0) {
      keep[i] <- TRUE
    }
  }
  if (drop) {
    if (any(!keep)) {
      fmt <- "Drop the following assay(s) for no data remained: %s"
      warn.msg <- sprintf(fmt, paste(which(drop), collapse = ", "))
      warning(warn.msg, immediate. = TRUE, call. = TRUE)
    }
    assays <- assays[keep]
    x@rowMap <- rowmap[, keep, drop = FALSE]
    x@colMap <- colmap[, keep, drop = FALSE]
  } else {
    x@rowMap <- rowmap
    x@colMap <- colmap
  }
  setRawAssays(x, assays, check = FALSE)
})

