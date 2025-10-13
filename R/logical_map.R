#' @include utils.R
NULL

#' @importFrom SparseArray cbind
#' @export
SparseArray_cbind <- SparseArray::cbind

#' @importFrom SparseArray rbind
#' @export
SparseArray_rbind <- SparseArray::rbind

#' Create a Sparse Logical Map
#'
#' Create a sparse logical matrix indicating the presence of each value across
#' a set of observations. Each row corresponds to a \code{\link{character}}
#' value, stored in row names. Each column represents an observation.
#'
#' @param x A character vector. Values to map across observations.
#' @param ncol Number of column (observation) for the output matrix.
#' @param names Names of observations. Must have `ncol` elements.
#' @param ... `r .dot_param`
#'
#' @returns
#' An `r .doc_links("SVT_SparseMatrix")` with `length(x)` rows and and `ncol`
#' columns. The matrix entries are logical, indicating whether each value in
#' `x` is present in each observation.
#'
#' @details
#' This helper is typically used internally within \code{\link{FlexAssays}} or
#' similar containers to track shared or mapped features across multiple
#' datasets or samples, while minimizing memory usage via sparse representation.
#'
#' @examples
#' # Create a sparse logical map
#' map <- sparseLogMap(LETTERS)
#'
#' # Create a logical map for known observations
#' obs <- letters[1:4]
#' map <- sparseLogMap(LETTERS, ncol = length(obs), names = obs)
#' map
#'
#' @importFrom SparseArray SVT_SparseArray
#'
#' @rdname sparseLogMap
#' @export
sparseLogMap <- function(x = character(), ncol = 0, names = NULL, ...) {
  if (missing(ncol)) {
    ncol <- length(names)
  }
  nrow <- length(x)
  out <- SVT_SparseArray(dim = c(nrow, ncol), type = "logical")
  colnames(out) <- names
  rownames(out) <- x
  out
}

#' Helper Functions for Sparse Logical Maps
#'
#' A set of utility functions for manipulating sparse logical maps, typically
#' used to track row mappings across observations. Includes methods for
#' accessing or modifying mappings, selecting intersected features, and dropping
#' unused rows or columns.
#'
#' @param x A sparse logical map, typically a `r .doc_links("SVT_SparseMatrix")`
#' @param ... `r .dot_param`
#'
#' @examples
#' # Create a logical map for known observations
#' obs <- letters[1:4]
#' map <- sparseLogMap(LETTERS[1:10], ncol = length(obs), names = obs)
#' map
#'
#' @name logmap-helpers
NULL

#' @param i Index or name of the column to access or modify in `x`.
#' @param value A new mapping to assign to column `i`. Must be one of the
#' followings:
#' \itemize{
#' \item A character vector: entries in `value` are set to `TRUE` for column
#' `i`.
#' \item `NULL`: all entries in column `i` are set to `FALSE`.
#' }
#' @param append Logical. If `TRUE`, values not currently in `rownames(x)` will
#' be appended as new rows. If `FALSE`, adding non-existing values will raise
#' an error.
#'
#' @returns
#' \itemize{
#' \item `mappedRowNames<-`: An updated logical map with modified mappings of
#' column `i`.
#' }
#'
#' @examples
#' # Set entries to a new column
#' mappedRowNames(map, 5) <- sample(LETTERS, 10)
#' mappedRowNames(map, "e") <- sample(LETTERS, 10)
#'
#' # Set entries to an existing column
#' mappedRowNames(map, "c") <- sample(LETTERS, 10)
#'
#' # Set 'c' to NULL. All values in column 'c' will be FALSE
#' mappedRowNames(map, "c") <- NULL
#' map
#'
#' @rdname logmap-helpers
#' @export
"mappedRowNames<-" <- function(x, i, append = TRUE, ..., value) {
  validLogMap(x)
  i <- i[1]
  if (is.character(i)) {
    idx <- match(i, colnames(x))
    if (is.na(idx)) {
      idx <- ncol(x) + 1L
    }
    return(.add_mapped_row_names(x, idx, value, new.name = i, append = append))
  }
  i <- as.integer(i)
  .add_mapped_row_names(x, i, value, append = append, ...)
}

.add_mapped_row_names <- function(x, i, value, new.name = NULL, append = TRUE) {
  if (i > ncol(x) + 1) {
    stop("Subscript out of bounds: ", i)
  }
  new.mat <- logical(nrow(x))
  not.found <- NULL
  if (length(value) > 0) {
    if (!is.character(value)) {
      stop("New row names must be a character vector.")
    }
    if (anyDuplicated(value)) {
      stop("New row names must not be duplicated.")
    }
    idx <- match(value, rownames(x), nomatch = 0L)
    new.mat[idx] <- TRUE
    not.found <- value[idx == 0]
  }
  if (i <= ncol(x)) {
    x[, i] <- new.mat
  } else {
    new.mat <- SVT_SparseArray(new.mat, dim = c(nrow(x), 1), type = "logical")
    colnames(new.mat) <- new.name
    if (ncol(x) > 0) {
      x <- SparseArray_cbind(x, new.mat)
    } else {
      rownames(new.mat) <- rownames(x)
      x <- new.mat
    }
  }
  if (length(not.found) == 0) {
    return(x)
  }
  if (!append) {
    stop("Cannot add non-existing rownames for when 'append' is FALSE")
  }
  empty.mat <- sparseLogMap(not.found, ncol(x), colnames(x))
  empty.mat[, i] <- TRUE
  SparseArray_rbind(x, empty.mat)
}

#' @param invert Logical. If `TRUE`, returns row names that are **not** mapped
#' to column `i`.
#'
#' @returns
#' \itemize{
#' \item `mappedRowNames`: A character vector specifying row names of `x` mapped
#' to column `i`.
#' }
#'
#' @examples
#' # Get row names that exist in a specific column
#' mappedRowNames(map, 1)
#' mappedRowNames(map, "e", invert = TRUE)
#'
#' @rdname logmap-helpers
#' @export
mappedRowNames <- function(x, i, invert = FALSE, ...) {
  validLogMap(x)
  if (invert) {
    return(rownames(x)[which(!x[, i, drop = TRUE])])
  }
  rownames(x)[which(x[, i, drop = TRUE])]
}

#' @param type Type of row intersection to return:
#' \itemize{
#' \item `all`(the default): features shared by all columns.
#' \item `duplicated`: features shared by at least **two** columns.
#' }
#'
#' @returns
#' \itemize{
#' \item `intersectedRows`: An integer vector of row indices mapped to multiple
#' columns, depending on `type`.
#' }
#'
#' @examples
#' # Which rows contain names that are mapped to multiple columns.
#' intersectedRows(map)
#'
#' @importFrom SparseArray rowSums
#'
#' @rdname logmap-helpers
#' @export
intersectedRows <- function(x, type = c("all", "duplicated"), ...) {
  validLogMap(x)
  type <- match.arg(type)
  if (type == "all") {
    return(which(rowSums(x) == ncol(x)))
  }
  which(rowSums(x) > 1)
}

#' @returns
#' \itemize{
#' \item `dropMapRows`: An updated map `x` with rows that are not mapped to any
#' column removed.
#' }
#'
#' @examples
#' # Drop rows whose names have no column mapped to.
#' dropMapRows(map)
#'
#' @importFrom SparseArray rowSums
#'
#' @rdname logmap-helpers
#' @export
dropMapRows <- function(x, ...) {
  validLogMap(x)
  fidx <- which(rowSums(x) == 0)
  if (length(fidx) > 0) {
    x <- x[-fidx, , drop = FALSE]
  }
  x
}

#' @returns
#' \itemize{
#' \item `dropMapCols`: An updated map `x` with columns that do not map to any
#' row removed.
#' }
#'
#' @examples
#' # Drop columns that no row names are mapped to.
#' dropMapCols(map)
#'
#' @importFrom SparseArray colSums
#'
#' @rdname logmap-helpers
#' @export
dropMapCols <- function(x, ...) {
  validLogMap(x)
  fidx <- which(colSums(x) == 0)
  if (length(fidx) > 0) {
    x <- x[, -fidx, drop = FALSE]
  }
  x
}

#' @returns
#' \itemize{
#' \item `removeMapCols`: An updated `x` with column `i` removed.
#' }
#'
#' @examples
#' # Remove columns
#' removeMapCols(map, 1)
#' removeMapCols(map, c(1, 3))
#' removeMapCols(map, c("a", "e"))
#'
#' @rdname logmap-helpers
#' @export removeMapCols
removeMapCols <- function(x, i, ...) {
  validLogMap(x)
  if (is.character(i)) {
    .valid_logical_map_colnames(x)
    old_i <- i
    i <- match(i, colnames(x))
    if (anyNA(i)) {
      stop("Subscript out of bounds: ", paste(old_i, collapse = ", "))
    }
  }
  .remove_map_cols(x, as.integer(i))
}

.remove_map_cols <- function(x, i) {
  if (ncol(x) == 0) {
    return(x)
  }
  i <- i[i > 0]
  if (length(i) == 0) {
    return(x)
  }
  if (any(i > ncol(x))) {
    stop("Subscript out of bounds (", ncol(x), "): ", i[i > ncol(x)][[1]])
  }
  x[, -i, drop = FALSE]
}

#' @param immediate. Logical. If `TRUE` (default), validation errors are thrown
#' immediately. If `FALSE`, errors are collected silently and return. This can
#' be useful in context of `r .doc_links("validObject")`.
#'
#' @returns
#' \itemize{
#' \item `validLogMap`: Check whether `x` is a valid sparse logical map.
#' }
#'
#' @rdname logmap-helpers
#' @export
validLogMap <- function(x, immediate. = TRUE) {
  invisible(c(
    .valid_classes_logical_map(x, immediate.),
    .valid_logical_map_colnames(x, immediate.),
    .valid_logical_map_rownames(x, immediate.)
  ))
}

#' @importClassesFrom SparseArray SVT_SparseMatrix
.valid_classes_logical_map <- function(x, immediate. = TRUE) {
  if (inherits(x, c("SVT_SparseMatrix"))) {
    return(invisible(NULL))
  }
  fmt <- "Invalid logical map class: <%s>."
  getErrors(sprintf(fmt, class(x)[1]), immediate. = immediate.)
}

.valid_logical_map_colnames <- function(x, immediate. = TRUE) {
  if (ncol(x) == 0) {
    return(invisible(NULL))
  }
  if (length(colnames(x)) == 0) {
    return(invisible(NULL))
  }
  if (!anyDuplicated(colnames(x), incomparables = "")) {
    return(invisible(NULL))
  }
  getErrors("Duplicate column names not allowed.", immediate. = immediate.)
}

.valid_logical_map_rownames <- function(x, immediate. = TRUE) {
  if (nrow(x) == 0) {
    return(invisible(NULL))
  }
  # Check rownames
  rnames <- rownames(x)
  err <- NULL
  if (length(rnames) == 0) {
    e <- "Row names must be supplied for a logical map."
    err <- getErrors(e, err, immediate.)
  }
  if (any(!nzchar(rnames))) {
    e <- "Row names cannot be empty strings for a logical map."
    err <- getErrors(e, err, immediate.)
  }
  if (anyDuplicated(rnames)) {
    e <- "Duplicate row names not allowed for a logical map."
    err <- getErrors(e, err, immediate.)
  }
  invisible(err)
}
