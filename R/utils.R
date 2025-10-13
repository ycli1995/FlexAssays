
#' Collect and handle errors in a standardized way
#'
#' `getErrors()` helps accumulate error messages during validity checks
#' or other internal validations. It supports returning the accumulated
#' errors or stopping immediately.
#'
#' @param e A character vector of new error messages to report.
#' @param err A character vector of existing errors (can be `NULL` or empty).
#' @param immediate. Logical. If `TRUE` (default), throws an error immediately
#' with the contents of `e`. If `FALSE`, returns all errors invisibly.
#'
#' @return If `immediate.` is `FALSE`, returns a combined character vector of
#' `err` and `e` invisibly. If `immediate.` is `TRUE`, raises an error.
#'
#' @examples
#' # Accumulate errors without stopping
#' err <- getErrors("First error", immediate. = FALSE)
#' err <- getErrors("Second error", err = err, immediate. = FALSE)
#' print(err)
#'
#' # getErrors("Something went wrong")
#'
#' @export
getErrors <- function(e = NULL, err = NULL, immediate. = TRUE) {
  if (length(e) == 0) {
    return(invisible(err))
  }
  if (immediate.) {
    stop(e, call. = FALSE)
  }
  invisible(c(err, e))
}

#' Display the S4 class title for an object
#'
#' This is a helper function to print the class name of an S4 object in
#' a standard format. It is usually called within \code{\link{show}} methods.
#'
#' @param object An S4 object.
#' @param ... Additional strings to include in the output.
#'
#' @return No return value. This function is called for printing.
#'
#' @examples
#' setClass("NewClass", slots = c(data = "character"))
#' showS4Title(new("NewClass"), "with character")
#'
#' @export
showS4Title <- function(object, ...) {
  cat(sprintf("%s object", class(object)[1]), ..., "\n")
}

#' Compare Two Objects Ignoring Specific Attributes
#'
#' This function compares two R objects using \code{\link{identical}}, but
#' allows users to ignore specific attributes during comparison. It is
#' especially useful when structural identity matters more than meta data.
#'
#' @param x First object to compare.
#' @param y Second object to compare.
#' @param ignore.attrs A character vector of attribute names to ignore during
#' comparison.
#' @param ... Additional arguments passed to \code{\link{identical}}.
#'
#' @return A logical scalar. `TRUE` when `x` and `y` are identical (ignoring
#' the specified attributes).
#'
#' @examples
#' a <- matrix(1:4, 2)
#' b <- matrix(1:4, 2)
#' attr(a, "meta") <- "something"
#' identical(a, b)    # FALSE
#' identicalNoAttrs(a, b, "meta")   # TRUE
#'
#' @export
identicalNoAttrs <- function(x, y, ignore.attrs = NULL, ...) {
  for (attr in ignore.attrs) {
    attr(x, attr) <- NULL
    attr(y, attr) <- NULL
  }
  identical(x = x, y = y, ...)
}

.identical_fmatch <- function(x, y, ...) {
  identicalNoAttrs(x, y, ignore.attrs = ".match.hash", ...)
}

#' Assign Dimnames with Sanity Checks
#'
#' A helper function for safely assigning dimnames to a matrix-like object.
#' Primarily used in method dispatch (e.g., `dimnames<-` for `FlexAssays`).
#'
#' @param x A matrix-like object with two dimensions.
#' @param value A list of length 2, containing the new row and column names
#' respectively.
#'
#' @return The modified object \code{x} with updated dimnames.
#'
#' @details
#' If the input \code{value} is not a list of length 2, will raise an error.
#' If the list is empty but the object \code{x} has non-zero rows or columns,
#' a warning is issued and the original object is returned unmodified to
#' prevents silent loss of dimension names on non-empty objects.
#'
#' This function is intended to be used inside `setMethod("dimnames<-", ...)`
#' and not called directly by end-users.
#'
#' @examples
#' mat <- matrix(1:4, nrow = 2)
#' autoDimnames(mat, list(c("A", "B"), c("X", "Y")))
#'
#' # Remove all dimension names
#' autoDimnames(mat, list(NULL, NULL))
#'
#' @export
autoDimnames <- function(x, value) {
  if (!(is.list(value) && length(value) == 2L)) {
    stop("dimnames replacement value must be a list of length 2")
  }
  rownames(x) <- value[[1]]
  colnames(x) <- value[[2]]
  x
}

#' Fast Intersection Between Two Vectors
#'
#' Performs a fast intersection between two vectors, with optional control
#' over whether to preserve duplicated elements in the first input.
#'
#' @param x,y Vectors (of the same mode) containing a sequence of items.
#' @param keep.duplicated Logical. If `TRUE`, duplicates in `x` are preserved.
#' Otherwise, duplicates are removed before matching.
#'
#' @return A vector of elements from `x` that are also in `y`. The result
#' preserves the order in `y`. Returns `NULL` if either input is `NULL`.
#'
#' @examples
#' fastIntersect(c("a", "b", "b", "c"), c("b", "c", "d"))
#' # [1] "b" "b" "c"
#'
#' fastIntersect(c("b", "b", "c"), c("b", "c", "d"), keep.duplicated = FALSE)
#' # [1] "b" "c"
#'
#' fastIntersect(NULL, c("a", "b"))
#' # NULL
#'
#' @export
fastIntersect <- function(x, y, keep.duplicated = TRUE) {
  if (is.null(x) || is.null(y)) {
    return(NULL)
  }
  if (keep.duplicated) {
    x <- unique(x)
  }
  ind <- match(x, y, nomatch = 0L)
  y[ind]
}

#' Reset Row and Column Names of a Matrix-like Object
#'
#' Replaces the row and/or column names of a matrix or matrix-like object.
#'
#' @param mat A matrix-like object. Must support \code{\link{rownames}} and
#' \code{\link{colnames}} methods.
#' @param rownames A character vector for new row names, or `NULL` to remove
#' row names.
#' @param colnames A character vector for new column names, or `NULL` to remove
#' column names.
#'
#' @return The input matrix with updated row and/or column names. The default
#' arguments will clear original row names and column names.
#'
#' @examples
#' mat <- matrix(1:4, nrow = 2)
#' resetDimNames(mat, rownames = c("g1", "g2"), colnames = c("c1", "c2"))
#' resetDimNames(mat) # returns mat with no row/column names
#
#' @export
resetDimNames <- function(mat, rownames = NULL, colnames = NULL) {
  rownames(mat) <- rownames
  colnames(mat) <- colnames
  mat
}

#' Subset a Matrix by Row and Column Names
#'
#' Efficiently subsets a matrix or matrix-like object by provided row and column
#' names. If the input names are already in the same order as the matrix, the
#' subsetting will be ignored to improve performance.
#'
#' @param mat A matrix-like object with row and column names.
#' @param rownames A character vector of row names to subset. Use zero-length
#' character to remove all rows.
#' @param colnames A character vector of column names to subset. Use zero-length
#' character to remove all columns.
#' @param drop Logical. Whether or not to drop dimensions (passed to matrix
#' subsetting).
#'
#' @return A subset of `mat` restricted to the provided `rownames` and
#' `colnames`.
#'
#' @details
#' If the input row/column names exactly match those in `mat`, the function
#' avoids lookup and directly returns the original matrix or slices it
#' efficiently.
#'
#' @examples
#' mat <- matrix(1:9, nrow = 3)
#' rownames(mat) <- c("geneA", "geneB", "geneC")
#' colnames(mat) <- c("cell1", "cell2", "cell3")
#'
#' subsetMatByDimNames(mat, rownames = "geneA", colnames = c("cell2", "cell3"))
#' # Returns a 1x2 matrix with values from geneA and specified columns
#'
#' @export
subsetMatByDimNames <- function(mat, rownames, colnames, drop = TRUE) {
  check.rows <- length(rownames) > 0
  check.cols <- length(colnames) > 0
  if (!check.rows) {
    if (!check.cols) {
      return(mat[integer(), integer(), drop = drop])
    }
    return(mat[integer(), colnames, drop = drop])
  }
  if (!check.cols) {
    return(mat[rownames, integer(), drop = drop])
  }
  check.rows <- .identical_fmatch(rownames(mat), rownames)
  check.cols <- .identical_fmatch(colnames(mat), colnames)
  if (check.rows) {
    if (!check.cols) {
      mat <- mat[, colnames, drop = drop]
    }
    return(mat)
  }
  if (check.cols) {
    return(mat[rownames, , drop = drop])
  }
  mat[rownames, colnames, drop = drop]
}
