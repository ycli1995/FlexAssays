#' @include utils.R
#' @include logical_map.R
#' @include FlexAssays-development.R
#' @importFrom methods as callNextMethod is new setClass
#' @importFrom utils head
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class
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
#' @importFrom methods setClass
#'
#' @importClassesFrom methods VIRTUAL
#' @importClassesFrom S4Vectors Annotated
#' @importClassesFrom SparseArray SVT_SparseMatrix
#'
#' @name FlexAssays
#' @docType class
#' @aliases SimpleFlexAssays-class FlexAssays-class
#'
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

#' @importFrom methods callNextMethod
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
#' @importClassesFrom S4Vectors SimpleList
#'
#' @exportClass SimpleFlexAssays
#' @rdname FlexAssays
setClass(
  "SimpleFlexAssays",
  contains = "FlexAssays",
  slots = c(data = "SimpleList")
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Constructor
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
#' @importFrom S4Vectors SimpleList
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
  .new_SimpleFlexAssays(
    assays,
    rnames = rownames,
    cnames = colnames,
    rowFlex = rowFlex,
    colFlex = colFlex,
    assayClasses = assayClasses,
    ...
  )
}

#' @importFrom S4Vectors new2
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

#' @importFrom S4Vectors new2
.new_SimpleFlexAssays <- function(
    assays,
    rnames = NULL,
    cnames = NULL,
    rowFlex = "free",
    colFlex = "free",
    assayClasses = NULL,
    ...
) {
  if (missing(assays)) {
    assays <- NULL
  }
  assays <- .normarg_assays(assays)
  if (length(assays) == 0) {
    return(.empty_SimpleFlexAssays(
      rnames = rnames,
      cnames = cnames,
      rowFlex = rowFlex,
      colFlex = colFlex,
      assayClasses = assayClasses
    ))
  }
  if (length(rnames) > 0 | length(cnames) > 0) {
    stop("'rownames' and 'colnames' only work when 'assays' is empty.")
  }

  .valid_assays_classes(assays, assayClasses)
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

#' @importFrom methods as
#' @importFrom S4Vectors new2
#' @importClassesFrom S4Vectors List SimpleList
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Validation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom S4Vectors setValidity2
setValidity2("FlexAssays", function(x) .validFlexAssays(x, immediate. = FALSE))

.validFlexAssays <- function(x, immediate. = TRUE) {
  err <- c(
    validLogMap(x@rowMap, immediate. = immediate.),
    validLogMap(x@colMap, immediate. = immediate.)
  )

  assays <- try(assays(x, withDimnames = FALSE), silent = TRUE)
  if (inherits(assays, "try-error")) {
    e <- "'assays(x)' must work for a FlexAssays."
    return(getErrors(e, err, immediate. = immediate.))
  }
  if (!is(assays, "SimpleList")) {
    e <- "'assays(x)' must return a SimpleList object."
    return(getErrors(e, err, immediate. = immediate.))
  }
  if (length(assays) == 0L) {
    return(invisible(err))
  }
  err <- c(
    .valid_mapped_assays(assays, x@rowMap, x@colMap, immediate. = immediate.),
    .valid_assays_classes(assays, x@assayClasses, immediate. = immediate.)
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
      err <- getErrors(e, err, immediate. = immediate.)
    }
    if (ncol(assays[[i]]) != mapped.ncols[i]) {
      e <- sprintf(fmt, "ncol", i, nrow(assays[[i]]), mapped.nrows[i])
      err <- getErrors(e, err, immediate. = immediate.)
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
  invisible(err)
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
  if (!anyNA(match(dimfun(mat), ref))) {
    return(invisible(NULL))
  }
  fmt <- "Not allowed '%s(%s)' found for 'bounded' flex: %s..."
  bad.names <- paste(head(setdiff(dimfun(mat), ref), 3), collapse = ", ")
  e <- sprintf(fmt, .get_func_name(dimfun), mat.name, bad.names)
  getErrors(e, immediate. = immediate.)
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
    ref <- dimfun(assays[[1]])
  }
  err <- NULL
  for (i in seq_along(assays)) {
    mat.name <- paste0("assays[[", i, "]]")
    err <- c(err, .valid_dnames_fixed(
      mat = assays[[i]],
      dimfun = dimfun,
      ref = ref,
      mat.name = mat.name,
      immediate. = immediate.
    ))
  }
  invisible(err)
}

.valid_dnames_fixed <- function(
    mat,
    dimfun,
    ref,
    mat.name = "",
    immediate. = TRUE
) {
  if (setequal(ref, dimfun(mat))) {
    return(invisible(NULL))
  }
  funstr <- .get_func_name(dimfun)

  bad.names <- head(setdiff(ref, dimfun(mat)), 3)
  if (length(bad.names) > 0) {
    fmt <- "Missed '%s(%s)' for 'fixed' flex: %s..."
    e <- sprintf(fmt, funstr, mat.name, paste(bad.names, collapse = ", "))
    return(getErrors(e, immediate. = immediate.))
  }

  bad.names <- head(setdiff(dimfun(mat), ref), 3)
  fmt <- "Out-of-bound '%s(%s)' for 'fixed' flex: %s..."
  e <- sprintf(fmt, funstr, mat.name, paste(bad.names, collapse = ", "))
  getErrors(e, immediate. = immediate.)
}

.valid_mapped_assays <- function(assays, rowmap, colmap, immediate. = TRUE) {
  e <- .valid_assays_length(assays, rowmap, colmap, immediate. = immediate.)
  if (length(assays) == 0) {
    return(invisible(e))
  }
  err <- e
  e <- .valid_assays_names(assays, rowmap, colmap, immediate. = immediate.)
  invisible(c(err, e))
}

.valid_assays_length <- function(assays, rowmap, colmap, immediate. = TRUE) {
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

.valid_assays_names <- function(assays, rowmap, colmap, immediate. = TRUE) {
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

.get_func_name <- function(func) {
  as.character(substitute(func))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Methods for FlexAssays
#'
#' A collection of accessor methods for retrieving or modifying basic structural
#' information from a \code{\link{FlexAssays}} object.
#'
#' @param x,object A \code{\link{FlexAssays}} object
#' @param ... `r .dot_param`
#' @param value `r .val_param` Typically a character vector for names or
#' dimension names.
#'
#' @examples
#' # Create some sparse matrices
#' fa <- exampleFlexAssays()
#' fa
#'
#' @name FlexAssays-methods
NULL

#' @returns
#' \itemize{
#' \item `show`: Print the basic information of `object` to the console.
#' }
#'
#' @note
#' These methods are primarily intended for interacting with the top-level
#' structure of a \code{\link{FlexAssays}} object.
#'
#' @importFrom S4Vectors coolcat
#' @importFrom methods show
#'
#' @export
#' @rdname FlexAssays-methods
setMethod("show", "FlexAssays", function(object) {
  showS4Title(object, sprintf("%ix%i", nrow(object), ncol(object)))
  cat(.show_flex(object))
  cat("assayClasses: ", paste(object@assayClasses, collapse = ", "), "\n")
  coolcat("rownames (%d): %s\n", rownames(object))
  coolcat("colnames (%d): %s\n", colnames(object))
  cat(sprintf("assays (%i):\n", length(object)))
  if (length(object) == 0) {
    return(invisible(NULL))
  }
  assay.names <- names(object)
  if (length(assay.names) == 0) {
    assay.names <- as.character(seq_along(object))
  }
  assays <- assays(object, withDimnames = FALSE)
  for (i in seq_along(assays)) {
    cat(sprintf(
      "  %s: %i x %i %s\n",
      assay.names[i],
      nrow(assays[[i]]),
      ncol(assays[[i]]),
      class(assays[[i]])[1]
    ))
  }
  invisible(NULL)
})

.show_flex <- function(x) {
  sprintf("rowFlex: '%s' colFlex: '%s'\n", x@rowFlex, x@colFlex)
}

#' @returns
#' \itemize{
#' \item `length`: Returns an integer indicating the number of assays.
#' }
#'
#' @examples
#' length(fa)
#'
#' @export
#' @rdname FlexAssays-methods
setMethod("length", "FlexAssays", function(x) {
  length(assays(x, withDimnames = FALSE))
})

#' @returns
#' \itemize{
#' \item `names`: A character vector specifying the layer names.
#' }
#'
#' @examples
#' names(fa)
#'
#' @export
#' @rdname FlexAssays-methods
setMethod("names", "FlexAssays", function(x) {
  names(assays(x, withDimnames = FALSE))
})

#' @returns
#' \itemize{
#' \item `names<-`: Returns `x` with updated assay names.
#' }
#'
#' @examples
#' names(fa) <- paste0("assay", as.character(seq_along(fa)))
#'
#' @export
#' @rdname FlexAssays-methods
setMethod("names<-", "FlexAssays", function(x, value) {
  colnames(x@rowMap) <- colnames(x@colMap) <- value
  assays <- assays(x, withDimnames = FALSE)
  names(assays) <- value
  setRawAssays(x, assays, check = FALSE)
})

#' @returns
#' \itemize{
#' \item `dim`: Returns an integer vector of length 2 representing the global
#' dimensions (rows × columns).
#' }
#'
#' @examples
#' dim(fa)
#' nrow(fa)
#' ncol(fa)
#'
#' @export
#' @rdname FlexAssays-methods
setMethod("dim", "FlexAssays", function(x) c(nrow(x@rowMap), nrow(x@colMap)))

#' @returns
#' \itemize{
#' \item `dimnames`: Returns a list of global row and column names.
#' }
#'
#' @export
#' @rdname FlexAssays-methods
setMethod("dimnames", "FlexAssays", function(x) {
  list(rownames(x@rowMap), rownames(x@colMap))
})

#' @returns
#' \itemize{
#' \item `dimnames<-`, `rownames<-` and `colnames<-`: Return an updated `x` with
#' modified row/column names.
#' }
#'
#' @export
#' @rdname FlexAssays-methods
setMethod("dimnames<-", "FlexAssays", function(x, value) {
  if (length(value) == 0 & any(dim(x) > 0)) {
    fmt <- "dimnames(x) cannot be empty for <%s> with non-empty rows or colums."
    warning(sprintf(fmt, class(x)[1]), immediate. = TRUE, call. = FALSE)
    return(x)
  }
  autoDimnames(x, value)
})

#' @importFrom BiocGenerics rownames<-
#'
#' @export
#' @rdname FlexAssays-methods
setMethod("rownames<-", "FlexAssays", function(x, value) {
  if (length(value) == 0 & nrow(x) > 0) {
    fmt <- "Cannot set empty row names to <%s> with non-empty rows."
    warning(sprintf(fmt, class(x)[1]), immediate. = TRUE, call. = FALSE)
    return(x)
  }
  rownames(x@rowMap) <- value
  .valid_logical_map_rownames(x@rowMap)
  x
})

#' @importFrom BiocGenerics colnames<-
#'
#' @export
#' @rdname FlexAssays-methods
setMethod("colnames<-", "FlexAssays", function(x, value) {
  if (length(value) == 0 & ncol(x) > 0) {
    fmt <- "Cannot set empty column names to <%s> with non-empty columns."
    warning(sprintf(fmt, class(x)[1]), immediate. = TRUE, call. = FALSE)
    return(x)
  }
  rownames(x@colMap) <- value
  .valid_logical_map_rownames(x@colMap)
  x
})

#' @returns
#' \itemize{
#' \item `as.list`: Returns a list of assays extracted from `x`.
#' }
#'
#' @rdname FlexAssays-methods
#' @export
#' @method as.list FlexAssays
as.list.FlexAssays <- function(x, ...) as.list(assays(x, ...))

#' @importFrom methods coerce setAs
setAs("SimpleList", "SimpleFlexAssays", function(from) {
  .new_SimpleFlexAssays(from)
})

setAs("SimpleList", "FlexAssays", function(from) {
  .new_SimpleFlexAssays(from)
})

setAs("FlexAssays", "SimpleList", function(from) {
  assays(from, withDimnames = TRUE)
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# assays
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get or set assays
#'
#' Methods to access or update the assay list for a  \code{\link{FlexAssays}}.
#'
#' @param x A \code{\link{FlexAssays}} object.
#' @param ... `r .dot_param`
#'
#' @examples
#' ## assays
#' fa <- exampleFlexAssays()
#'
#' @name FlexAssays-assays
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
#' @details
#' \itemize{
#' \item `assays` method must be implemented for each concrete sub-class of
#' `FlexAssays`.
#' }
#'
#' @examples
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
})

#' @param value An object of a class specified in the S4 method signature.
#' \itemize{
#' \item `assays<-`: A list of matrix-like objects to set as the internal assay
#' data of `x`. If is `NULL`, all assays in `x` will be removed.
#' \item `assay<-`: A single matrix-like object to set as the `i`-th assay. Must
#' contain both row and column names. If is `NULL`, the original assay at index
#' `i` will be removed.
#' }
#'
#' @details
#' [`assays<-`] method is defined for the VIRTUAL class [`FlexAssays`] and
#' should automatically ensure the output `x` is validated. Any sub-class that
#' implements [`setRawAssays`] and [`cleanRawAssays`] correctly should works for
#' [`assays<-`].
#'
#' @seealso [`setRawAssays`] and [`cleanRawAssays`]
#'
#' @examples
#' ## assays<-
#' fa <- exampleFlexAssays()
#' aa <- assays(fa)[1:2]
#' assays(fa) <- aa
#'
#' @rdname FlexAssays-assays
#' @export
setMethod(
  "assays<-", c(x = "FlexAssays", value = "list_OR_List"),
  function(x, withDimnames = TRUE, ..., value) {
    x <- cleanRawAssays(x)
    x <- .set_assays(x, seq_along(value), value)
    names(x) <- names(value)
    x
  }
)

#' @examples
#' ## Remove all assays
#' fa <- exampleFlexAssays()
#' assays(fa) <- NULL
#'
#' @rdname FlexAssays-assays
#' @export
setMethod(
  "assays<-", c(x = "FlexAssays", value = "NULL"),
  function(x, withDimnames = TRUE, ..., value) {
    cleanRawAssays(x)
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
      fmt <- "assay(<%s>, i=\"missing\") failed:\n  No assay exists. "
      stop(sprintf(fmt, class(x)[1]))
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
  fmt <- "'assay(<%s>, i=\"%s\", ...)' invalid subscript 'i': %s\n"
  msg <- sprintf(fmt, class(x)[1], class(i)[1], as.character(i))
  assay <- tryCatch(
    assays(x, withDimnames = FALSE)[[i]],
    error = function(err) stop(msg, conditionMessage(err))
  )
  if (is.null(assay)) {
    stop(msg, sprintf("Non-existing assay '%s'", as.character(i)))
  }
  if (!withDimnames) {
    return(assay)
  }
  resetDimNames(
    mat = assay,
    rownames = mappedRowNames(x@rowMap, i = i),
    colnames = mappedRowNames(x@colMap, i = i)
  )
}

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

.set_assays <- function(x, which, new.assays) {
  assays <- assays(x, withDimnames = FALSE)
  if (is.null(new.assays)) {  # Remove some assays (or all of them)
    x@rowMap <- removeMapCols(x@rowMap, which)
    x@colMap <- removeMapCols(x@colMap, which)
    assays[which] <- NULL
    return(setRawAssays(x, assays, check = FALSE))
  }
  if (is.data.frame(new.assays)) {
    new.assays <- list(new.assays)
  }
  if (!inherits(new.assays, c("list", "List"))) {
    new.assays <- list(new.assays)
  }
  stopifnot(length(which) == length(new.assays))
  for (ii in seq_along(which)) {  # Add or update selected assay
    i <- which[ii]
    new.assay <- new.assays[[ii]]
    .valid_new_assay(x = x, i = i, new.assay = new.assay)
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
  }
  setRawAssays(x, assays, check = FALSE)
}

.set_one_assay <- function(x, i, new.assay) {
  assays <- assays(x, withDimnames = FALSE)
  if (is.null(new.assay)) {
    x@rowMap <- removeMapCols(x@rowMap, i)
    x@colMap <- removeMapCols(x@colMap, i)
    assays[[i]] <- new.assay
    return(setRawAssays(x, assays, check = FALSE))
  }
  .valid_new_assay(x = x, i = i, new.assay = new.assay)
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
  setRawAssays(x, assays, check = FALSE)
}

.valid_new_assay <- function(x, i, new.assay) {
  .valid_assay_classes(new.assay, x@assayClasses)
  if (nrow(new.assay) > 0 & length(rownames(new.assay)) == 0) {
    stop("New assay '", i, "' must have row names with non-empty rows.")
  }
  if (ncol(new.assay) > 0 & length(colnames(new.assay)) == 0) {
    stop("New assay '", i, "' must have column names with non-empty columns.")
  }
  if (x@rowFlex == "fixed") {
    .valid_dnames_fixed(
      new.assay,
      dimfun = rownames,
      ref = rownames(x@rowMap),
      mat.name = paste0("<assays[[", i, "]]>")
    )
  }
  if (x@rowFlex == "bounded") {
    .valid_dnames_bounded(
      new.assay,
      dimfun = rownames,
      ref = rownames(x@rowMap),
      mat.name = paste0("<assays[[", i, "]]>")
    )
  }
  if (x@colFlex == "fixed") {
    .valid_dnames_fixed(
      new.assay,
      dimfun = colnames,
      ref = rownames(x@colMap),
      mat.name = paste0("<assays[[", i, "]]>")
    )
  }
  if (x@colFlex == "bounded") {
    .valid_dnames_bounded(
      new.assay,
      dimfun = colnames,
      ref = rownames(x@colMap),
      mat.name = paste0("<assays[[", i, "]]>")
    )
  }
  invisible(NULL)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Example data
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
  FlexAssays(list(m1 = m1, m2, m3 = m3))
}

.generate_sparse_matrix <- function(nrow, ncol, density = 0.75) {
  mat1 <- Matrix::rsparsematrix(nrow = nrow, ncol = ncol, density = density)
  rownames(mat1) <- sample(LETTERS, nrow(mat1))
  colnames(mat1) <- sample(letters, ncol(mat1))
  mat1
}
