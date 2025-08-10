#' @include logical_map.R
#' @include utils.R
#'
#' @importFrom S4Vectors coolcat isTRUEorFALSE List new2 setValidity2 SimpleList
#' @importFrom Matrix rsparsematrix
#' @importFrom methods as callNextMethod coerce is new setAs setClass
#' @importFrom utils .DollarNames head
#' @importClassesFrom methods VIRTUAL
#' @importClassesFrom S4Vectors character_OR_NULL SimpleList
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class ########################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The FlexAssays class
#'
#' The `FlexAssays` virtual class provide a formal abstraction of layered
#' matrix-like objects (`assays`) with potential common names of dimensions.
#' Rows and columns of each assay are mapped using sparse logical maps (`rowMap`
#' and `colMap`), so that users can track whether an entry exists in rows or
#' columns for a specific assay.
#'
#' @slot rowMap,colMap An `r .doc_links("SVT_SparseMatrix")` to map row (column)
#' names for each assay.
#'
#' @slot rowFlex,colFlex See the description for parameters in the constructor
#' function.
#'
#' @slot assayClasses See the description for parameters in the constructor
#' function.
#'
#' @details
#' For a given `FlexAssays`, `rowFlex`, `colFlex` and `assayClasses` must only
#' be set when the instance is initialized. The end-users should not modify
#' these slots.
#'
#' @name FlexAssays
#' @docType class
#' @aliases SimpleFlexAssays-class FlexAssays-class
#' @exportClass FlexAssays
setClass(
  Class = "FlexAssays",
  contains = c("Annotated", "VIRTUAL"),
  slots = c(
    assayClasses = "character_OR_NULL",
    rowMap = "SVT_SparseMatrix",
    colMap = "SVT_SparseMatrix",
    rowFlex = "character",
    colFlex = "character"
  )
)

setMethod("initialize", "FlexAssays", function(
    .Object, ...,
    rowFlex = c("free", "bounded", "fixed"),
    colFlex = c("free", "bounded", "fixed")
) {
  .Object@rowFlex <- match.arg(rowFlex)
  .Object@colFlex <- match.arg(colFlex)
  callNextMethod(.Object, ...)
})

#' @slot data In `SimpleFlexAssays`, `data` represents a list storing all
#' matrix-like assays. Row names of each assay must exist in row names of
#' `rowMap`, and so do column names for `colMap`.
#'
#' @exportClass SimpleFlexAssays
#' @rdname FlexAssays
setClass(
  "SimpleFlexAssays",
  contains = "FlexAssays",
  slots = c(data = "SimpleList")
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Constructor ##################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param assays One of the following:
#' \itemize{
#' \item A single matrix-like object
#' \item A list of multiple matrix-like objects.
#' \item An existing `FlexAssays`.
#' }
#' @param ... `r .dot_param`
#' @param rownames,colnames A character vector containing the expected row or
#' column names. Default is `NULL`, which will fetch the union of row/column
#' names in `assays`. Only used when `assays` is empty or an existing
#' `FlexAssays`.
#' @param assayClasses A character vector indicating the supported classes of
#' assays. Default is `NULL` indicating no restrict of assay class.
#' @param rowFlex,colFlex A character scalar indicating the flexibility of rows
#' or columns. Must be one of the following:
#' \itemize{
#' \item `free`: The global rows/columns can be extended when a new assay with
#' non-existing dimension names is added.
#' \item `bounded`: The global rows/columns **cannot** be extended when a new
#' assay with non-existing dimension names is added.
#' \item `fixed`: All assays must contain the same row/column names.
#' }
#'
#' @seealso
#' [rowFlex()] and [colFlex()] to learn how to control the flexibility of
#' dimensions.
#'
#' @examples
#' # Create some sparse matrices
#' mat1 <- Matrix::rsparsematrix(nrow = 10, ncol = 15, density = 0.75)
#' dimnames(mat1) <- list(sample(LETTERS, 10), sample(letters, 15))
#'
#' mat2 <- Matrix::rsparsematrix(nrow = 15, ncol = 12, density = 0.75)
#' dimnames(mat2) <- list(sample(LETTERS, 15), sample(letters, 12))
#'
#' mat3 <- Matrix::rsparsematrix(nrow = 25, ncol = 10, density = 0.75)
#' dimnames(mat3) <- list(sample(LETTERS, 25), sample(letters, 10))
#'
#' # Create a FlexAssays from a list
#' mats <- FlexAssays(list(m1 = mat1, m2 = mat2, m3 = mat3))
#' mats
#'
#' # Create a FlexAssays from a single matrix
#' mats <- FlexAssays(mat1)
#' mats
#'
#' # Create a FlexAssays from an existing FlexAssays
#' mats2 <- FlexAssays(
#'   mats,
#'   rownames = head(rownames(mats)),
#'   colnames = head(colnames(mats))
#' )
#'
#' @rdname FlexAssays
#' @export FlexAssays
FlexAssays <- function(
    assays = SimpleList(),
    rownames = NULL,
    colnames = NULL,
    assayClasses = NULL,
    rowFlex = "free",
    colFlex = "free",
    ...
) {
  if (inherits(assays, "FlexAssays")) {
    return(assays)
  }
  assays <- .normarg_assays(assays)
  .init_SimpleFlexAssays(
    assays,
    rnames = rownames,
    cnames = colnames,
    rowFlex = rowFlex,
    colFlex = colFlex,
    assayClasses = assayClasses,
    ...
  )
}

.normarg_assays <- function(assays) {
  err <- paste(
    "'assays' must be a list or SimpleList of matrix-like elements,",
    "or a matrix-like object, or a NULL."
  )
  if (is.null(assays)) {
    return(SimpleList())
  }
  if (length(dim(assays)) == 2L) {
    return(new2("SimpleList", listData = list(assays), check = FALSE))
  }
  if (!is(assays, "SimpleList")) {
    if (is.list(assays)) {
      return(new2("SimpleList", listData = assays, check = FALSE))
    }
    if (is(assays, "List")) {
      return(as(assays, "SimpleList"))
    }
    stop(err)
  }
  assays
}

.empty_SimpleFlexAssays <- function(
    rnames = NULL,
    cnames = NULL,
    rowFlex = "free",
    colFlex = "free",
    assayClasses = character()
) {
  new2(
    Class = "SimpleFlexAssays",
    check = FALSE,
    rowMap = sparseLogMap(as.character(rnames)),
    colMap = sparseLogMap(as.character(cnames)),
    rowFlex = rowFlex,
    colFlex = colFlex,
    assayClasses = assayClasses
  )
}

.init_SimpleFlexAssays <- function(
    assays,
    rnames = NULL,
    cnames = NULL,
    rowFlex = "free",
    colFlex = "free",
    assayClasses = NULL,
    ...
) {
  if (length(assays) == 0) {
    return(.empty_SimpleFlexAssays(
      rnames = rnames,
      cnames = cnames,
      rowFlex = rowFlex,
      colFlex = colFlex,
      assayClasses = assayClasses
    ))
  }
  .valid_assays_classes(assays, assayClasses)
  if (length(rnames) > 0 | length(cnames) > 0) {
    stop("'rownames' and 'colnames' only work when 'assays' is empty.")
  }
  if (rowFlex == "fixed") {
    .valid_assays_dnames_fixed(assays, dimfun = rownames)
  }
  if (colFlex == "fixed") {
    .valid_assays_dnames_fixed(assays, dimfun = colnames)
  }
  rnames <- unique(unlist(lapply(assays, rownames)))
  cnames <- unique(unlist(lapply(assays, colnames)))
  rowMap <- sparseLogMap(as.character(rnames), length(assays), names(assays))
  colMap <- sparseLogMap(as.character(cnames), length(assays), names(assays))
  for (i in seq_along(assays)) {
    mappedRowNames(rowMap, i) <- rownames(assays[[i]])
    mappedRowNames(colMap, i) <- colnames(assays[[i]])
    assays[[i]] <- resetDimNames(subsetMatByDimNames(
      mat = assays[[i]],
      rownames = mappedRowNames(rowMap, i),
      colnames = mappedRowNames(colMap, i),
      drop = FALSE
    ))
  }
  new2(
    Class = "SimpleFlexAssays",
    check = FALSE,
    data = assays,
    rowMap = rowMap,
    colMap = colMap,
    rowFlex = rowFlex,
    colFlex = colFlex,
    assayClasses = assayClasses,
    ...
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Validation ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setValidity2("FlexAssays", function(x) validFlexAssays(x, immediate. = FALSE))

#' Validation for FlexAssays
#'
#' Check whether a \code{\link{FlexAssays}} is valid.
#'
#' @param x A \code{\link{FlexAssays}} to check.
#' @param immediate. Logical indicating if the errors should be raised
#' immediately. Set `immediate. = FALSE` can be useful for context of
#' `r .doc_links("validObject")`.
#'
#' @seealso [getErrors()]
#'
#' @returns
#' Returns an invisible `NULL` if `x` is a valid \code{\link{FlexAssays}}. If
#' `immediate. = TRUE`, this function will raise an error when `x` is invalid.
#' Otherwise, it will return a character vector that collects all errors.
#'
#' @examples
#' m1 <- matrix(1:10, 2, 5)
#' rownames(m1) <- letters[1:2]
#' colnames(m1) <- LETTERS[1:5]
#'
#' fa <- FlexAssays(m1)
#' validFlexAssays(fa)
#'
#' @export
validFlexAssays <- function(x, immediate. = TRUE) {
  err <- c(validLogMap(x@rowMap, immediate.), validLogMap(x@colMap, immediate.))

  assays <- try(assays(x, withDimnames = FALSE), silent = TRUE)
  if (inherits(assays, "try-error")) {
    e <- "'assays(x)' must work for a FlexAssays."
    return(getErrors(e, err, immediate.))
  }
  if (!is(assays, "SimpleList")) {
    e <- "'assays(x)' must return a SimpleList object."
    return(getErrors(e, err, immediate.))
  }
  if (length(assays) == 0L) {
    return(invisible(err))
  }
  err <- c(
    .valid_mapped_assays_list(assays, x@rowMap, x@colMap, immediate.),
    .valid_assays_classes(assays, x@assayClasses, immediate.)
  )
  ## Check dimension numbers for each assay.
  if (x@rowFlex == "fixed") {
    mapped.nrows <- rep.int(nrow(x@rowMap), ncol(x@rowMap))
  } else {
    mapped.nrows <- colSums(x@rowMap)
  }
  if (x@colFlex == "fixed") {
    mapped.ncols <- rep.int(nrow(x@colMap), ncol(x@colMap))
  } else {
    mapped.ncols <- colSums(x@colMap)
  }
  fmt <- paste0(.show_flex(x), " `%s(assays[[%d]])` (%d) is not equal to %d`")
  for (i in seq_along(assays)) {
    if (nrow(assays[[i]]) != mapped.nrows[i]) {
      e <- sprintf(fmt, "nrow", i, nrow(assays[[i]]), mapped.nrows[i])
      err <- getErrors(e, err, immediate.)
    }
    if (ncol(assays[[i]]) != mapped.ncols[i]) {
      e <- sprintf(fmt, "ncol", i, nrow(assays[[i]]), mapped.nrows[i])
      err <- getErrors(e, err, immediate.)
    }
  }
  invisible(err)
}

.valid_assay_classes <- function(assay, classes = NULL, immediate. = TRUE) {
  if (length(classes) == 0) {
    return(invisible(NULL))
  }
  if (inherits(assay, classes)) {
    return(invisible(NULL))
  }
  fmt <- "Invalid assay <%s>. Valid 'assayClasses' include: %s"
  e <- sprintf(fmt, class(assay)[1], paste(classes, collapse = ", "))
  getErrors(e, immediate. = immediate.)
}

.valid_assays_classes <- function(assays, classes = NULL, immediate. = TRUE) {
  if (length(classes) == 0) {
    return(invisible(NULL))
  }
  err <- NULL
  for (i in seq_along(assays)) {
    e <- .valid_assay_classes(assays[[i]], classes, immediate. = immediate.)
    err <- c(err, e)
  }
  return(invisible(err))
}

.valid_dnames_bounded <- function(
    mat,
    dimfun,
    ref = NULL,
    mat.name = "",
    immediate. = TRUE
) {
  if (length(ref) == 0) {
    return(invisible(NULL))
  }
  if (anyNA(match(dimfun(mat), ref))) {
    fmt <- "Not allowed '%s(%s)' found for 'bounded' flex: %s..."
    funstr <- as.character(substitute(dimfun))
    bad.names <- paste(head(setdiff(dimfun(mat), ref), 3), collapse = ", ")
    return(getErrors(
      sprintf(fmt, funstr, mat.name, bad.names),
      immediate. = immediate.
    ))
  }
  return(invisible(NULL))
}

.valid_dnames_fixed <- function(
    mat,
    dimfun,
    ref = NULL,
    mat.name = "",
    immediate. = TRUE
) {
  if (length(ref) == 0) {
    return(invisible(NULL))
  }
  if (setequal(ref, dimfun(mat))) {
    return(invisible(NULL))
  }
  funstr <- as.character(substitute(dimfun))
  bad.names <- head(setdiff(ref, dimfun(mat)), 3)
  if (length(bad.names) > 0) {
    fmt <- "Missed '%s(%s)' for 'fixed' flex: %s..."
    return(getErrors(
      sprintf(fmt, funstr, mat.name, paste(bad.names, collapse = ", ")),
      immediate. = immediate.
    ))
  }
  bad.names <- paste(head(setdiff(dimfun(mat), ref), 3), collapse = ", ")
  fmt <- "Out-of-bound '%s(%s)' for 'fixed' flex: %s..."
  getErrors(sprintf(fmt, funstr, mat.name, bad.names), immediate. = immediate.)
}

.valid_assays_dnames_fixed <- function(
    assays,
    dimfun,
    ref = NULL,
    immediate. = TRUE
) {
  if (length(ref) == 0) {
    if (length(assays) == 1) {
      return(invisible(NULL))
    }
    if (Reduce(setequal, lapply(assays, dimfun))) {
      return(invisible(NULL))
    }
    funstr <- as.character(substitute(dimfun))
    fmt <- "All assays must contain the same '%s()'"
    return(getErrors(sprintf(fmt, funstr), immediate. = immediate.))
  }
  err <- NULL
  for (i in seq_along(assays)) {
    err <- c(err, .valid_dnames_fixed(
      mat = assays[[i]],
      dimfun = dimfun,
      ref = ref,
      mat.name = paste0("assays[[", i, "]]"),
      immediate. = immediate.
    ))
  }
  invisible(err)
}

.valid_assays_dnames_bounded <- function(
    assays,
    dimfun,
    ref = NULL,
    immediate. = TRUE
) {
  err <- NULL
  for (i in seq_along(assays)) {
    err <- c(err, .valid_dnames_bounded(
      mat = assays[[i]],
      dimfun = dimfun,
      ref = ref,
      mat.name = paste0("assays[[", i, "]]"),
      immediate. = immediate.
    ))
  }
  invisible(err)
}

.valid_flex_assay_one_dim <- function(
    new.assay, dimmap, flex,
    by = c("row", "col"),
    immediate. = TRUE,
    ...
) {
  if (flex == "free") {
    return(invisible(NULL))
  }
  by <- match.arg(by)
  if (flex == "fixed") {
    return(.valid_dnames_fixed(
      mat = new.assay,
      dimfun = if (by == "row") rownames else colnames,
      ref = rownames(dimmap),
      immediate. = immediate.,
      ...
    ))
  }
  .valid_dnames_bounded(
    mat = new.assay,
    dimfun = if (by == "row") rownames else colnames,
    ref = rownames(dimmap),
    immediate. = immediate.,
    ...
  )
}

.valid_flex_assays_one_dim <- function(
    new.assays, dimmap, flex,
    by = c("row", "col"),
    immediate. = TRUE
) {
  if (flex == "free") {
    return(invisible(NULL))
  }
  by <- match.arg(by)
  if (flex == "fixed") {
    return(.valid_assays_dnames_fixed(
      assays = new.assays,
      dimfun = if (by == "row") rownames else colnames,
      ref = rownames(dimmap),
      immediate. = immediate.
    ))
  }
  .valid_assays_dnames_bounded(
    assays = new.assays,
    dimfun = if (by == "row") rownames else colnames,
    ref = rownames(dimmap),
    immediate. = immediate.
  )
}

.valid_mapped_assays_length <- function(
    assays,
    rowmap,
    colmap,
    immediate. = TRUE
) {
  n.assays <- length(assays)
  err <- NULL
  if (ncol(rowmap) != n.assays) {
    e <- "`ncol(rowMap)` must be equal to `length(assays)`."
    err <- getErrors(e, err, immediate.)
  }
  if (ncol(colmap) != n.assays) {
    e <- "`ncol(colMap)` must be equal to `length(assays)`."
    err <- getErrors(e, err, immediate.)
  }
  invisible(err)
}

.valid_mapped_assays_names <- function(
    assays,
    rowmap,
    colmap,
    immediate. = TRUE
) {
  rowmap.names <- colnames(rowmap)
  colmap.names <- colnames(colmap)
  assay.names <- names(assays)
  err <- NULL
  if (length(assay.names) == 0) {
    fmt <- "Invalid `colnames(%s)` for empty assay name: %s"
    if (length(rowmap.names) > 0) {
      e <- sprintf(fmt, "rowMap", paste(rowmap.names, collapse = ", "))
      err <- getErrors(e, err, immediate.)
    }
    if (length(colmap.names) > 0) {
      e <- sprintf(fmt, "colMap", paste(colmap.names, collapse = ", "))
      err <- getErrors(e, err, immediate.)
    }
    return(invisible(err))
  }
  fmt <- "Unmatched assay names with `colnames(%s)`: %s"
  if (!all(assay.names == rowmap.names)) {
    e <- sprintf(fmt, "rowMap", paste(assay.names, collapse = ", "))
    err <- getErrors(e, err, immediate.)
  }
  if (!all(assay.names == colmap.names)) {
    e <- sprintf(fmt, "colMap", paste(assay.names, collapse = ", "))
    err <- getErrors(e, err, immediate.)
  }
  invisible(err)
}

.valid_mapped_assays_list <- function(
    assays,
    rowmap,
    colmap,
    immediate. = TRUE
) {
  err <- .valid_mapped_assays_length(assays, rowmap, colmap, immediate.)
  if (length(assays) == 0) {
    return(invisible(err))
  }
  err <- c(err, .valid_mapped_assays_names(assays, rowmap, colmap, immediate.))
  invisible(err)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Coerce #######################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setAs("SimpleList", "SimpleFlexAssays", function(from) {
  .init_SimpleFlexAssays(from)
})

setAs("SimpleList", "FlexAssays", function(from) {
  .init_SimpleFlexAssays(from)
})

setAs("FlexAssays", "SimpleList", function(from) {
  assays(from, withDimnames = TRUE)
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# example data #################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' A Mini Example for FlexAssays
#'
#' Generate a tiny \code{\link{FlexAssays}} object as example.
#'
#' @returns
#' A \code{\link{FlexAssays}} object with 3 assays.
#'
#' @examples
#' exampleFlexAssays()
#'
#' @export
exampleFlexAssays <- function() {
  m1 <- .generate_sparse_matrix(10, 6)
  m2 <- .generate_sparse_matrix(8, 13)
  m3 <- .generate_sparse_matrix(7, 5)
  FlexAssays(list(m1, m2, m3))
}

.generate_sparse_matrix <- function(nrow, ncol, density = 0.75) {
  mat1 <- rsparsematrix(nrow = nrow, ncol = ncol, density = density)
  rownames(mat1) <- sample(LETTERS, nrow(mat1))
  colnames(mat1) <- sample(letters, ncol(mat1))
  return(mat1)
}


