#' @include FlexAssays-class.R
#' @include FlexAssays-development.R
#'
#' @importFrom BiocGenerics colnames<- rownames<-
#' @importFrom S4Vectors coolcat isTRUEorFALSE List new2 setValidity2 SimpleList
#' @importFrom methods show
#' @importFrom utils .DollarNames
#' @importClassesFrom methods VIRTUAL
#' @importClassesFrom S4Vectors character_OR_NULL list_OR_List SimpleList
NULL

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
#' mat1 <- Matrix::rsparsematrix(nrow = 20, ncol = 15, density = 0.75)
#' dimnames(mat1) <- list(sample(LETTERS, 20), sample(letters, 15))
#'
#' mat2 <- Matrix::rsparsematrix(nrow = 15, ncol = 12, density = 0.75)
#' dimnames(mat2) <- list(sample(LETTERS, 15), sample(letters, 12))
#'
#' mat3 <- Matrix::rsparsematrix(nrow = 25, ncol = 10, density = 0.75)
#' dimnames(mat3) <- list(sample(LETTERS, 25), sample(letters, 10))
#'
#' mats <- FlexAssays(list(m1 = mat1, m2 = mat2, m3 = mat3))
#' mats
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
#' structure of a \code{\link{FlexAssays}} object. Some methods (e.g.,
#' `showFlexAssays`) are meant to assist developers during inspection or
#' debugging.
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
#' length(mats)
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
#' names(mats)
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
#' @export
#' @rdname FlexAssays-methods
setMethod("names<-", "FlexAssays", function(x, value) {
  colnames(x@rowMap) <- colnames(x@colMap) <- value
  assays <- assays(x, withDimnames = FALSE)
  names(assays) <- value
  setRawAssays(x, assays)
})

#' @returns
#' \itemize{
#' \item `dim`: Returns an integer vector of length 2 representing the global
#' dimensions (rows × columns).
#' }
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

#' @export
#' @rdname FlexAssays-methods
setMethod("colnames<-", "FlexAssays", function(x, value) {
  if (length(value) == 0 & nrow(x) > 0) {
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
